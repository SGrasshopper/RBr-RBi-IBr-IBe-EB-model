# RBr-RBi-IBr-IBe-EB-model
Model of Chlamydial development

CellState.py is needed to replace the stock version to add additional variables to the model

Rbr matures into Rbe based on percentchance curve from empirical data

Rbe divides > Rbe and IB

IB starts with hctA == 0, matures to EB(celltype == 5), based on HctB accumulation                 
      
ISSUES WITH CURRENT MODEL
Curently gene expression controls cell type and this is good.  
But would like to make gene expression controlling gene expression 
[Euo] should repress [HctA] and [CtcB] and hctA feedback should repress [Euo]
[ctcB] should promote [hctB] and [hctB] should repress everthing.
