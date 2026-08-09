#!/usr/bin/env python3
"""Build the canonical portable-report input for the ICDAR benchmark."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


DEFAULT_PREFIX = "tree-construction-comparison"
COMPONENT_AND_TOS = {
    "tos_max4c_min8c",
    "tos_self_dual",
    "max_tree_8c",
    "min_tree_8c",
}
RESIDUAL = {
    "residual_unrestricted_8c",
    "residual_saturated_8c",
}
CHART_LABELS = {
    "tos_max4c_min8c": "ToS M4/M8",
    "tos_self_dual": "ToS SD",
    "max_tree_8c": "Max 8c",
    "min_tree_8c": "Min 8c",
    "residual_unrestricted_8c": "Irrestrita 8c",
    "residual_saturated_8c": "Saturada 8c",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("results_dir", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--generated-at", required=True)
    parser.add_argument("--prefix", default=DEFAULT_PREFIX)
    return parser.parse_args()


def read_csv(path: Path) -> list[dict[str, object]]:
    numeric_fields = {
        "pixels",
        "megapixels",
        "mean_ms",
        "mean_ci95_low_ms",
        "mean_ci95_high_ms",
        "median_ms",
        "standard_deviation_ms",
        "mean_ms_per_megapixel",
        "mean_nodes",
        "median_repetition_relative_range_pct",
        "rank",
        "relative_to_fastest",
        "empirical_exponent",
        "time_ratio_1080p_over_480p",
        "pixel_ratio_1080p_over_480p",
    }
    with path.open(newline="", encoding="utf-8") as stream:
        rows: list[dict[str, object]] = []
        for source_row in csv.DictReader(stream):
            row: dict[str, object] = dict(source_row)
            for field in numeric_fields.intersection(row):
                row[field] = float(str(row[field]))
            rows.append(row)
        return rows


def source(
    identifier: str,
    label: str,
    path: str,
    description: str,
    generated_at: str,
    *,
    filters: list[str] | None = None,
    metric_definitions: list[str] | None = None,
    sql: str | None = None,
) -> dict[str, object]:
    return {
        "id": identifier,
        "label": label,
        "path": path,
        "query": {
            "description": description,
            "executed_at": generated_at,
            "language": "SQL" if sql else ("CSV" if path.endswith(".csv") else "JSON"),
            **({"engine": "DuckDB", "sql": sql} if sql else {}),
            "filters": filters or [],
            "metric_definitions": metric_definitions or [],
        },
    }


def main() -> None:
    args = parse_args()
    prefix = args.prefix
    post_cleanup = "post-cleanup" in prefix
    report_title = "Tempo de construção de árvores morfológicas no ICDAR"
    if post_cleanup:
        report_title += " — pós-limpeza"
    summary = read_csv(args.results_dir / f"{prefix}-summary.csv")
    scaling = read_csv(args.results_dir / f"{prefix}-scaling.csv")
    summary_by_key = {
        (str(row["resolution"]), str(row["algorithm"])): row
        for row in summary
    }

    def mean(resolution: str, algorithm: str) -> float:
        return float(summary_by_key[(resolution, algorithm)]["mean_ms"])

    def format_pt(value: float, digits: int = 2) -> str:
        rendered = f"{value:,.{digits}f}"
        return rendered.replace(",", "_").replace(".", ",").replace("_", ".")

    max_tree_1080 = mean("1080p", "max_tree_8c")
    min_tree_1080 = mean("1080p", "min_tree_8c")
    tos_self_dual_1080 = mean("1080p", "tos_self_dual")
    tos_max_min_1080 = mean("1080p", "tos_max4c_min8c")
    residual_unrestricted_1080 = mean("1080p", "residual_unrestricted_8c")
    residual_saturated_1080 = mean("1080p", "residual_saturated_8c")
    residual_ratio_1080 = residual_saturated_1080 / residual_unrestricted_1080
    saturation_overheads = [
        (
            mean(resolution, "residual_saturated_8c")
            / mean(resolution, "residual_unrestricted_8c")
            - 1.0
        )
        * 100.0
        for resolution in ("480p", "720p", "1080p")
    ]
    tos_advantages = [
        (
            mean(resolution, "tos_max4c_min8c")
            / mean(resolution, "tos_self_dual")
            - 1.0
        )
        * 100.0
        for resolution in ("480p", "720p", "1080p")
    ]
    scaling_time_ratios = [
        float(row["time_ratio_1080p_over_480p"])
        for row in scaling
    ]
    scaling_exponents = [float(row["empirical_exponent"]) for row in scaling]

    headline = {
        "max_tree_1080_ms": max_tree_1080,
        "max_tree_1080_vs_480": (
            mean("1080p", "max_tree_8c") / mean("480p", "max_tree_8c")
        ),
        "tos_self_dual_1080_ms": tos_self_dual_1080,
        "tos_self_dual_vs_fastest": float(
            summary_by_key[("1080p", "tos_self_dual")]["relative_to_fastest"]
        ),
        "residual_unrestricted_1080_ms": residual_unrestricted_1080,
        "residual_unrestricted_vs_fastest": float(
            summary_by_key[("1080p", "residual_unrestricted_8c")][
                "relative_to_fastest"
            ]
        ),
        "residual_saturated_1080_ms": residual_saturated_1080,
        "residual_saturated_vs_unrestricted": residual_ratio_1080,
    }

    common_summary_metrics = [
        "Tempo por imagem: mediana de cinco construções; agregado: média das medianas das dez imagens.",
        "IC 95%: bootstrap descritivo da média sobre as dez imagens, com 20.000 reamostragens.",
        "Tempo normalizado: média em milissegundos dividida pelo número de megapixels.",
        "Relativo ao mais rápido: média do método dividida pela menor média na mesma resolução.",
    ]
    sources = [
        source(
            "headline",
            "Indicadores principais em 1080p",
            f"benchmarks/results/{prefix}-summary.csv",
            "Seleção dos quatro tempos principais e de seus comparadores a partir da síntese validada.",
            args.generated_at,
            metric_definitions=common_summary_metrics,
            sql=(
                "WITH s AS (SELECT * FROM read_csv_auto('benchmarks/results/"
                f"{prefix}-summary.csv'))\n"
                "SELECT\n"
                "  max(CASE WHEN resolution='1080p' AND algorithm='max_tree_8c' THEN mean_ms END) AS max_tree_1080_ms,\n"
                "  max(CASE WHEN resolution='1080p' AND algorithm='max_tree_8c' THEN mean_ms END) / max(CASE WHEN resolution='480p' AND algorithm='max_tree_8c' THEN mean_ms END) AS max_tree_1080_vs_480,\n"
                "  max(CASE WHEN resolution='1080p' AND algorithm='tos_self_dual' THEN mean_ms END) AS tos_self_dual_1080_ms,\n"
                "  max(CASE WHEN resolution='1080p' AND algorithm='tos_self_dual' THEN relative_to_fastest END) AS tos_self_dual_vs_fastest,\n"
                "  max(CASE WHEN resolution='1080p' AND algorithm='residual_unrestricted_8c' THEN mean_ms END) AS residual_unrestricted_1080_ms,\n"
                "  max(CASE WHEN resolution='1080p' AND algorithm='residual_unrestricted_8c' THEN relative_to_fastest END) AS residual_unrestricted_vs_fastest,\n"
                "  max(CASE WHEN resolution='1080p' AND algorithm='residual_saturated_8c' THEN mean_ms END) AS residual_saturated_1080_ms,\n"
                "  max(CASE WHEN resolution='1080p' AND algorithm='residual_saturated_8c' THEN mean_ms END) / max(CASE WHEN resolution='1080p' AND algorithm='residual_unrestricted_8c' THEN mean_ms END) AS residual_saturated_vs_unrestricted\n"
                "FROM s;"
            ),
        ),
        source(
            "summary",
            "Síntese do tempo de construção",
            f"benchmarks/results/{prefix}-summary.csv",
            "Agregados por resolução e algoritmo derivados das medianas por imagem.",
            args.generated_at,
            filters=[
                "test_000.png a test_009.png",
                "resoluções 480p, 720p e 1080p",
                "cinco repetições por imagem e algoritmo",
            ],
            metric_definitions=common_summary_metrics,
            sql=(
                "SELECT * FROM read_csv_auto('benchmarks/results/"
                f"{prefix}-summary.csv') ORDER BY pixels, rank;"
            ),
        ),
        source(
            "scaling",
            "Escalabilidade empírica",
            f"benchmarks/results/{prefix}-scaling.csv",
            "Ajuste log-log da média do tempo contra o número de pixels nas três resoluções.",
            args.generated_at,
            metric_definitions=[
                "Expoente empírico p: inclinação de log(tempo médio) contra log(número de pixels).",
                "Razão 1080p/480p: média em 1080p dividida pela média em 480p.",
            ],
            sql=(
                "SELECT * FROM read_csv_auto('benchmarks/results/"
                f"{prefix}-scaling.csv') ORDER BY empirical_exponent DESC;"
            ),
        ),
        source(
            "per_image",
            "Medianas e dispersão por imagem",
            f"benchmarks/results/{prefix}-per-image.csv",
            "Mediana, mínimo, máximo e amplitude relativa das cinco repetições de cada grupo.",
            args.generated_at,
        ),
        source(
            "validation",
            "Auditoria de completude",
            f"benchmarks/results/{prefix}-validation.json",
            "Resultado das verificações de cobertura, reconstrução, determinismo e positividade.",
            args.generated_at,
        ),
        source(
            "raw",
            "Medições brutas",
            f"benchmarks/results/{prefix}-raw.csv",
            "Novecentas medições individuais produzidas pelo executável C++ em modo Release.",
            args.generated_at,
        ),
    ]

    artifact = {
        "surface": "report",
        "manifest": {
            "version": 1,
            "surface": "report",
            "title": report_title,
            "description": (
                "Comparação de seis construções em dez imagens, três resoluções e "
                "cinco repetições, com medição externa da chamada pública completa."
            ),
            "generatedAt": args.generated_at,
            "sources": sources,
            "cards": [
                {
                    "id": "fastest-1080",
                    "description": "Menor média em 1080p; max-tree e min-tree ficaram praticamente empatadas.",
                    "dataset": "headline_1080p",
                    "sourceId": "headline",
                    "metrics": [
                        {
                            "label": "Max-tree 1080p (ms)",
                            "field": "max_tree_1080_ms",
                            "format": "number",
                        },
                        {
                            "label": "1080p / 480p",
                            "field": "max_tree_1080_vs_480",
                            "format": "number",
                        },
                    ],
                },
                {
                    "id": "tos-1080",
                    "description": "A ToS self-dual foi a ToS mais rápida em média nas três resoluções.",
                    "dataset": "headline_1080p",
                    "sourceId": "headline",
                    "metrics": [
                        {
                            "label": "ToS self-dual 1080p (ms)",
                            "field": "tos_self_dual_1080_ms",
                            "format": "number",
                        },
                        {
                            "label": "× tempo mínimo",
                            "field": "tos_self_dual_vs_fastest",
                            "format": "number",
                        },
                    ],
                },
                {
                    "id": "unrestricted-1080",
                    "description": "A residual irrestrita inclui a construção inicial da max-tree e da min-tree.",
                    "dataset": "headline_1080p",
                    "sourceId": "headline",
                    "metrics": [
                        {
                            "label": "Residual irrestrita 1080p (ms)",
                            "field": "residual_unrestricted_1080_ms",
                            "format": "number",
                        },
                        {
                            "label": "× tempo mínimo",
                            "field": "residual_unrestricted_vs_fastest",
                            "format": "number",
                        },
                    ],
                },
                {
                    "id": "saturated-1080",
                    "description": (
                        "A certificação de saturação adicionou cerca de "
                        f"{format_pt((residual_ratio_1080 - 1.0) * 100.0, 0)}% "
                        "em 1080p."
                    ),
                    "dataset": "headline_1080p",
                    "sourceId": "headline",
                    "metrics": [
                        {
                            "label": "Residual saturada 1080p (ms)",
                            "field": "residual_saturated_1080_ms",
                            "format": "number",
                        },
                        {
                            "label": "× irrestrita",
                            "field": "residual_saturated_vs_unrestricted",
                            "format": "number",
                        },
                    ],
                },
            ],
            "charts": [
                {
                    "id": "component-tos-times",
                    "title": "Tempo de construção das árvores de componentes e ToS",
                    "subtitle": "Média de dez medianas por imagem; cinco repetições por mediana; unidade: ms.",
                    "intent": "comparison",
                    "question": "Como max-tree, min-tree e as duas ToS se comparam nas três resoluções?",
                    "rationale": "Barras agrupadas preservam as três resoluções discretas e a comparação direta entre quatro métodos.",
                    "comparisonContext": {
                        "grain": "média de dez medianas por imagem",
                        "denominator": "dez imagens ICDAR",
                        "normalization": "cinco repetições por imagem e método",
                        "unit": "ms",
                    },
                    "type": "bar",
                    "dataset": "component_tos_summary",
                    "sourceId": "summary",
                    "encodings": {
                        "x": {
                            "field": "resolution",
                            "type": "ordinal",
                            "label": "Resolução",
                        },
                        "y": {
                            "field": "mean_ms",
                            "type": "quantitative",
                            "format": "number",
                            "label": "Tempo médio",
                            "unit": "ms",
                        },
                        "color": {
                            "field": "chart_label",
                            "type": "nominal",
                            "label": "Árvore",
                        },
                        "tooltip": [
                            {"field": "algorithm_label", "label": "Árvore"},
                            {"field": "mean_ci95_low_ms", "label": "IC 95% inferior"},
                            {"field": "mean_ci95_high_ms", "label": "IC 95% superior"},
                            {"field": "mean_ms_per_megapixel", "label": "ms/MP"},
                            {"field": "mean_nodes", "label": "Nós médios"},
                        ],
                    },
                    "valueFormat": "number",
                    "unit": "ms",
                    "layout": "full",
                    "labels": {"values": "all"},
                    "legend": {"position": "bottom", "title": "Árvore"},
                    "palette": {"kind": "categorical"},
                    "settings": {
                        "groupMode": "grouped",
                        "showValues": True,
                        "sort": "none",
                    },
                    "surface": {"surface": "card", "viewMode": "both"},
                    "maxRows": 12,
                },
                {
                    "id": "residual-times",
                    "title": "Tempo de construção das árvores residuais",
                    "subtitle": "Média de dez medianas por imagem; a versão saturada certifica a conectividade do complemento.",
                    "intent": "comparison",
                    "question": "Qual é o custo adicional da restrição de saturação nas três resoluções?",
                    "rationale": "Separar as residuais evita que sua escala comprima visualmente as árvores mais rápidas.",
                    "comparisonContext": {
                        "baseline": "residual irrestrita",
                        "grain": "média de dez medianas por imagem",
                        "denominator": "dez imagens ICDAR",
                        "normalization": "cinco repetições por imagem e método",
                        "unit": "ms",
                    },
                    "type": "bar",
                    "dataset": "residual_summary",
                    "sourceId": "summary",
                    "encodings": {
                        "x": {
                            "field": "resolution",
                            "type": "ordinal",
                            "label": "Resolução",
                        },
                        "y": {
                            "field": "mean_ms",
                            "type": "quantitative",
                            "format": "number",
                            "label": "Tempo médio",
                            "unit": "ms",
                        },
                        "color": {
                            "field": "chart_label",
                            "type": "nominal",
                            "label": "Árvore residual",
                        },
                        "tooltip": [
                            {"field": "algorithm_label", "label": "Árvore"},
                            {"field": "mean_ci95_low_ms", "label": "IC 95% inferior"},
                            {"field": "mean_ci95_high_ms", "label": "IC 95% superior"},
                            {"field": "mean_ms_per_megapixel", "label": "ms/MP"},
                            {"field": "mean_nodes", "label": "Nós médios"},
                        ],
                    },
                    "valueFormat": "number",
                    "unit": "ms",
                    "layout": "full",
                    "labels": {"values": "all"},
                    "legend": {"position": "bottom", "title": "Árvore residual"},
                    "palette": {"kind": "categorical"},
                    "settings": {
                        "groupMode": "grouped",
                        "showValues": True,
                        "sort": "none",
                    },
                    "surface": {"surface": "card", "viewMode": "both"},
                    "maxRows": 6,
                },
            ],
            "tables": [
                {
                    "id": "complete-results",
                    "title": "Resultados agregados completos",
                    "subtitle": "Média e IC 95% sobre dez medianas por imagem; tempo normalizado por megapixel.",
                    "dataset": "all_summary",
                    "defaultSort": {"field": "pixels", "direction": "asc"},
                    "density": "dense",
                    "sourceId": "summary",
                    "layout": "full",
                    "columns": [
                        {"field": "pixels", "label": "Pixels", "format": "compact"},
                        {"field": "resolution", "label": "Resolução", "type": "text"},
                        {"field": "algorithm_label", "label": "Árvore", "type": "text"},
                        {"field": "mean_ms", "label": "Média (ms)", "format": "number"},
                        {"field": "mean_ci95_low_ms", "label": "IC inf. (ms)", "format": "number"},
                        {"field": "mean_ci95_high_ms", "label": "IC sup. (ms)", "format": "number"},
                        {"field": "mean_ms_per_megapixel", "label": "ms/MP", "format": "number"},
                        {"field": "relative_to_fastest", "label": "× tempo mínimo", "format": "number"},
                    ],
                },
                {
                    "id": "scaling-results",
                    "title": "Escalabilidade entre 480p e 1080p",
                    "subtitle": "Expoente do ajuste log-log e razão observada de tempo; razão de pixels = 5,064×.",
                    "dataset": "scaling",
                    "defaultSort": {"field": "empirical_exponent", "direction": "desc"},
                    "density": "spacious",
                    "sourceId": "scaling",
                    "layout": "full",
                    "columns": [
                        {"field": "algorithm_label", "label": "Árvore", "type": "text"},
                        {"field": "empirical_exponent", "label": "Expoente p", "format": "number"},
                        {"field": "time_ratio_1080p_over_480p", "label": "Tempo 1080p/480p", "format": "number"},
                        {"field": "pixel_ratio_1080p_over_480p", "label": "Pixels 1080p/480p", "format": "number"},
                    ],
                },
            ],
            "blocks": [
                {
                    "id": "title",
                    "type": "markdown",
                    "layout": "full",
                    "body": f"# {report_title}",
                },
                {
                    "id": "technical-summary",
                    "type": "markdown",
                    "layout": "full",
                    "sourceId": "summary",
                    "body": (
                        "## Resumo técnico\n\n"
                        "As max-tree e min-tree 8c são as construções mais rápidas e ficam praticamente empatadas. "
                        f"Em 1080p, suas médias foram **{format_pt(max_tree_1080)} ms** e "
                        f"**{format_pt(min_tree_1080)} ms**. A ToS self-dual exigiu "
                        f"**{format_pt(tos_self_dual_1080)} ms** e a ToS Max4cMin8c, "
                        f"**{format_pt(tos_max_min_1080)} ms**. A residual irrestrita chegou a "
                        f"**{format_pt(residual_unrestricted_1080)} ms**; a saturada, a "
                        f"**{format_pt(residual_saturated_1080)} ms**, ou "
                        f"**{format_pt(residual_ratio_1080)}×** o tempo da irrestrita."
                    ),
                },
                {
                    "id": "headline-strip",
                    "type": "metric-strip",
                    "layout": "full",
                    "cardIds": [
                        "fastest-1080",
                        "tos-1080",
                        "unrestricted-1080",
                        "saturated-1080",
                    ],
                },
                {
                    "id": "component-finding",
                    "type": "markdown",
                    "layout": "full",
                    "sourceId": "summary",
                    "body": (
                        "## Árvores de componentes dominam; as duas ToS são próximas\n\n"
                        f"A max-tree e a min-tree custaram **{format_pt(max_tree_1080 / 1000.0)} s** e "
                        f"**{format_pt(min_tree_1080 / 1000.0)} s** em 1080p. As ToS custaram "
                        f"**{format_pt(min(tos_self_dual_1080, tos_max_min_1080) / 1000.0)}–"
                        f"{format_pt(max(tos_self_dual_1080, tos_max_min_1080) / 1000.0)} s**. "
                        "A diferença entre as duas ToS é pequena diante da diferença entre famílias; "
                        f"a vantagem relativa da mais rápida variou de **{format_pt(min(tos_advantages), 1)}%** "
                        f"a **{format_pt(max(tos_advantages), 1)}%** entre as resoluções."
                    ),
                },
                {
                    "id": "component-chart",
                    "type": "chart",
                    "layout": "full",
                    "chartId": "component-tos-times",
                    "sourceId": "summary",
                },
                {
                    "id": "residual-finding",
                    "type": "markdown",
                    "layout": "full",
                    "sourceId": "summary",
                    "body": (
                        "## A certificação de saturação é o maior sobrecusto\n\n"
                        f"A residual irrestrita foi **{format_pt(residual_unrestricted_1080 / tos_self_dual_1080)}×** "
                        "mais lenta que a ToS self-dual em 1080p. A residual saturada adicionou "
                        f"**{format_pt(min(saturation_overheads), 0)}–{format_pt(max(saturation_overheads), 0)}%** "
                        "sobre a irrestrita entre as três resoluções e, em 1080p, atingiu "
                        f"**{format_pt(residual_saturated_1080 / max_tree_1080)}×** o tempo da max-tree. "
                        "A dispersão entre imagens é muito maior nas residuais, "
                        "o que mostra dependência forte da estrutura de extremos e não apenas do número de pixels."
                    ),
                },
                {
                    "id": "residual-chart",
                    "type": "chart",
                    "layout": "full",
                    "chartId": "residual-times",
                    "sourceId": "summary",
                },
                {
                    "id": "complete-results-heading",
                    "type": "markdown",
                    "layout": "full",
                    "body": (
                        "## Valores completos confirmam a mesma ordenação\n\n"
                        "A tabela preserva médias, incerteza descritiva, custo por megapixel e razão contra o "
                        "método mais rápido em cada resolução."
                    ),
                },
                {
                    "id": "complete-results-table",
                    "type": "table",
                    "layout": "full",
                    "tableId": "complete-results",
                    "sourceId": "summary",
                },
                {
                    "id": "scaling-finding",
                    "type": "markdown",
                    "layout": "full",
                    "sourceId": "scaling",
                    "body": (
                        "## O crescimento observado é aproximadamente linear no número de pixels\n\n"
                        "De 480p para 1080p, o domínio cresce **5,064×** e os tempos crescem entre "
                        f"**{format_pt(min(scaling_time_ratios))}×** e **{format_pt(max(scaling_time_ratios))}×**. "
                        f"Os expoentes empíricos ficam entre **{format_pt(min(scaling_exponents))}** e "
                        f"**{format_pt(max(scaling_exponents))}**. Isso não é uma "
                        "prova de complexidade assintótica, mas indica que, neste intervalo, as diferenças são "
                        "principalmente fatores constantes específicos de cada construção."
                    ),
                },
                {
                    "id": "scaling-table",
                    "type": "table",
                    "layout": "full",
                    "tableId": "scaling-results",
                    "sourceId": "scaling",
                },
                {
                    "id": "scope",
                    "type": "markdown",
                    "layout": "full",
                    "sourceId": "raw",
                    "body": (
                        "## Escopo, dados e definições\n\n"
                        "Foram usadas `test_000.png`–`test_009.png` do ICDAR em **480×853**, **720×1280** "
                        "e **1080×1920**, sempre em tons de cinza de 8 bits. `Max4cMin8c` corresponde ao modo "
                        "interno `Min8cMax4c`, que nomeia primeiro a conectividade dos mínimos. Max-tree, "
                        "min-tree e as duas residuais usam adjacência 8c (`radius=1.5`); a residual saturada usa "
                        "o pixel de índice zero como infinito."
                    ),
                },
                {
                    "id": "methodology",
                    "type": "markdown",
                    "layout": "full",
                    "sourceId": "validation",
                    "body": (
                        "## Metodologia reproduzível\n\n"
                        "Cada método foi aquecido uma vez e validado por reconstrução exata. Depois foram feitas "
                        "cinco construções por imagem, com ordem rotacionada entre os seis métodos. O tempo mede "
                        "externamente a chamada pública completa e exclui carregamento da imagem e destruição da "
                        "árvore; o código de produção não coleta durações internas. Nas "
                        "residuais, inclui as max-tree e min-tree iniciais e a materialização final. O agregado "
                        "é a média das dez medianas por imagem. A execução foi Release/NDEBUG, sequencial, em "
                        "Apple M4 com Apple Clang 21."
                    ),
                },
                {
                    "id": "robustness",
                    "type": "markdown",
                    "layout": "full",
                    "sourceId": "per_image",
                    "body": (
                        "## Limitações, incerteza e robustez\n\n"
                        "A cobertura está completa: **900 medições**, sem tempo inválido, com número de nós "
                        "determinístico e reconstrução exata no aquecimento. Todas as repetições válidas foram "
                        "preservadas, e a mediana de cinco reduz a influência de atrasos isolados. Os ICs "
                        "bootstrap descrevem apenas a variação entre estas dez "
                        "imagens, que são um subconjunto fixo, não uma amostra aleatória de todo o ICDAR. Os "
                        "valores também são específicos desta máquina, compilador e versão do código; memória e "
                        "paralelismo não foram medidos. Como todas as famílias variaram entre as rodadas, a "
                        "comparação com séries anteriores não isola causalmente o efeito da limpeza de "
                        "instrumentação."
                    ),
                },
                {
                    "id": "recommendations",
                    "type": "markdown",
                    "layout": "full",
                    "body": (
                        "## Próximos passos recomendados\n\n"
                        "1. Priorizar o perfil da **certificação de saturação**, pois ela explica o maior "
                        "sobrecusto controlável, usando um profiler externo ou builds dedicadas.\n"
                        "2. Usar `test_003` e outras imagens lentas como casos de perfil, mantendo as dez imagens "
                        "para confirmar generalidade.\n"
                        "3. Preservar max-tree/min-tree 8c como piso de custo e repetir este protocolo após cada "
                        "otimização.\n"
                        "4. Acrescentar pico de memória e métricas de eventos apenas ao executável de benchmark "
                        "ou às implementações de referência, sem reintroduzir instrumentação na produção."
                    ),
                },
                {
                    "id": "further-questions",
                    "type": "markdown",
                    "layout": "full",
                    "body": (
                        "## Questões em aberto\n\n"
                        "Quanto da variação entre imagens é explicado pelo número de nós, extremos rejeitados e "
                        "travessias exatas do complemento? A ordem relativa permanece igual em imagens coloridas ou com "
                        "maior profundidade de níveis? E a equivalência de custo entre max-tree e min-tree se "
                        "mantém em outras distribuições de contraste?"
                    ),
                },
            ],
        },
        "snapshot": {
            "version": 1,
            "generatedAt": args.generated_at,
            "status": "ready",
            "datasets": {
                "headline_1080p": [headline],
                "component_tos_summary": [
                    {**row, "chart_label": CHART_LABELS[str(row["algorithm"])]}
                    for row in summary
                    if str(row["algorithm"]) in COMPONENT_AND_TOS
                ],
                "residual_summary": [
                    {**row, "chart_label": CHART_LABELS[str(row["algorithm"])]}
                    for row in summary
                    if str(row["algorithm"]) in RESIDUAL
                ],
                "all_summary": summary,
                "scaling": scaling,
            },
            "accessIssues": [],
        },
        "sources": sources,
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as stream:
        json.dump(artifact, stream, ensure_ascii=False, indent=2)
        stream.write("\n")


if __name__ == "__main__":
    main()
