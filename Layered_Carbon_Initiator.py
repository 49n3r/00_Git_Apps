import streamlit as st
import random as ran
import numpy as np
import pandas as pd
import time
import plotly
import plotly.express as px

### This Application builds random initial models for molecular dynamics simulation of 
#amorphous graphie, multi-shell fullerenes and carbon nanotubes

st.set_page_config(
    page_title="MTG-OU",
    page_icon="https://chinonsougwumadu.com/wp-content/uploads/2024/05/microsoftteams-image-17.jpg")

st.write("## Welcome to the MTG-OU App Center 👋")

st.sidebar.success("Select a constructor above.")
logo_url = "https://chinonsougwumadu.com/wp-content/uploads/2024/05/microsoftteams-image-17.jpg"
st.sidebar.image(logo_url)

st.markdown(
    """
    This app designed by the Materials Theory Group - Ohio University (MTG-OU) builds random 
    initial models for molecular dynamics (MD) simulation of amorphous graphite, multi-shell fullerenes
    and carbon nanotubes.

    **👈 Select a constructor from the sidebar** to build your initial models for MD simulations!
    ### Want to learn more aboout MTG-OU?
    - Check out our [resarch](https://daviddrabold.com/research/)
    - Jump into our [publications](https://daviddrabold.com/publications/)
    - Enjoy some cool animations on our [YouTube](https://youtube.com/playlist?list=PLkk5oFUpZDlbSOc7lLYdTramEX_Y3jx8m&si=UIAPj_2G2cCTeaUH) channel.
    ### Check out other MTG-OU apps
    - This [streamlit app](https://porous-carbon-constructor.streamlit.app/) (also included here 👈) constructs atomistic porous carbon models which 
      can be converted to continuum models for finite element simulation. See the associated publication here.
"""
)

