# %%
import celltypist
from celltypist import models
# %%
models.download_models(model = ['Immune_All_Low.pkl', 'Immune_All_High.pkl'], force_update = True)

# %%
model = celltypist.Model.load('Immune_All_Low.pkl')
model.convert()
model.write('Immune_All_Low_Mouse.pkl')
# %%
model = celltypist.Model.load('Immune_All_High.pkl')
model.convert()
model.write('Immune_All_High_Mouse.pkl')