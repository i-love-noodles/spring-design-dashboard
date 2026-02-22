import streamlit as st

st.set_page_config(page_title="Spring Design Analyzer", layout="wide")

pg = st.navigation([
    st.Page("pages/analysis.py", title="Analysis", default=True),
    st.Page("pages/optimization.py", title="Optimization"),
])
pg.run()
