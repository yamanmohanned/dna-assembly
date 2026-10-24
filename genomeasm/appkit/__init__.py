"""Everything the Streamlit application needs that is not drawing widgets.

The application (``app/``) is a thin layer: pages call these functions, which are tested without a
browser. Nothing here imports Streamlit. Charts are matplotlib figures so that every chart can be
exported as SVG and PNG with the same code, offline.
"""
