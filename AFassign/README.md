AFassign is a function which allows the user to assign anatomical features to each spot in the spatial dataset.

Input paramters:

obj: Seurat object of spatial data

method: the source tissue for spatial data i.e. FFPE/Frozen

thres: Threshold correlation value for AF assignment. The default is 0.15

<b>Example</b>

<i>AFassign(obj, 'FFPE', 0.15)</i>
