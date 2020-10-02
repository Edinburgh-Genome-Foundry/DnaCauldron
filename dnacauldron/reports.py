from datetime import datetime
import os

from pdf_reports import (
    write_report,
    pug_to_html,
)

from .version import __version__

THIS_PATH = os.path.dirname(os.path.realpath(__file__))
ASSETS_PATH = os.path.join(THIS_PATH, "report_assets")
DOMESTICATION_REPORT_TEMPLATE = os.path.join(ASSETS_PATH, "domestication_report.pug")
STYLESHEET = os.path.join(ASSETS_PATH, "report_style.css")


def dnacauldron_pug_to_html(template, **context):
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


def write_pdf_domestication_report(target):
    html = dnacauldron_pug_to_html(DOMESTICATION_REPORT_TEMPLATE,)
    write_report(html, target, extra_stylesheets=(STYLESHEET,))
