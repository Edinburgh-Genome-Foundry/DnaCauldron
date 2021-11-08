from datetime import datetime
import os

try:
    import pdf_reports

    PDF_REPORTS_AVAILABLE = True
except ImportError:
    PDF_REPORTS_AVAILABLE = False
else:
    from pdf_reports import (
        dataframe_to_html,
        style_table_rows,
        add_css_class,
        pug_to_html,
        write_report,
    )

from .version import __version__

THIS_PATH = os.path.dirname(os.path.realpath(__file__))
ASSETS_PATH = os.path.join(THIS_PATH, "report_assets")
DOMESTICATION_REPORT_TEMPLATE = os.path.join(ASSETS_PATH, "simulation_report.pug")
STYLESHEET = os.path.join(ASSETS_PATH, "report_style.css")


def dnacauldron_pug_to_html(template, **context):
    if not PDF_REPORTS_AVAILABLE:
        raise ImportError("Could not load PDF Reports")
    now = datetime.now().strftime("%Y-%m-%d")
    defaults = {
        "sidebar_text": "Generated on %s by DNA Cauldron version %s"
        % (now, __version__),
        "dc_logo_url": os.path.join(ASSETS_PATH, "imgs", "logo.png"),
    }
    for k in defaults:
        if k not in context:
            context[k] = defaults[k]
    return pug_to_html(template, **context)


def write_simulation_pdf_report(target, simulation_info):
    if not PDF_REPORTS_AVAILABLE:
        raise ImportError("Could not load PDF Reports")
    summary_table = dataframe_to_html(simulation_info, extra_classes=("definition",))

    def tr_modifier(tr):
        tds = list(tr.find_all("td"))
        if len(tds) == 0:
            return
        outcome, number = tds
        if outcome.text == "Valid":
            if number.text == "0":
                add_css_class(tr, "negative")
            else:
                add_css_class(tr, "positive")
        elif number.text != "0":
            add_css_class(tr, "negative")

    summary_table = style_table_rows(summary_table, tr_modifier)
    html = dnacauldron_pug_to_html(
        DOMESTICATION_REPORT_TEMPLATE, summary_table=summary_table
    )
    write_report(html, target, extra_stylesheets=(STYLESHEET,))
