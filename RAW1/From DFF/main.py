import pandas as pd

df_main = pd.read_excel("bees.xlsx", sheet_name="main")
df_site = pd.read_excel("bees.xlsx", sheet_name="nbr_sp_by_site.xls")

df3 = df_site.merge(
    df_main[["Species", "Statut_IUCN_Belgium"]],
    left_on="TAXPRIO",
    right_on="Species",
    how="left").drop(
    columns=["IUCN", "Species"])

df3 = df3.rename(columns={
    "Statut_IUCN_Belgium":"IUNC"
})

df3.to_excel("clean.xlsx", index=False)
