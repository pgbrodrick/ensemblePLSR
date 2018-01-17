# ensemblePLSR
The purpose of this code is to take a generically formatted csv database of chems and spectra, and run PLSR calibrations on that data. While multiple methods of running this analysis are viable using this code (through use of the associated settings file), this is adapted very directly from the second chapter of Dana Chadwick's thesis (Chadwick, K. Dana, and Gregory P. Asner. "Organismic-scale remote sensing of canopy foliar traits in lowland tropical forests." Remote Sensing 8.2 (2016): 87.). Chadwick's original code used the autopls package, which showed strong results, but is very slow to execute. This code adapts
large portions of that work in a more clearly documented, easier to modify, and faster way. The paper for the autopls package is (Schmidtlein, Sebastian, Hannes Feilhauer, and Helge Bruelheide. "Mapping plant strategy types using remote sensing." Journal of Vegetation Science 23.3 (2012): 395-405.).

More specifically, this code uses crown-based subsets of chemistry and spectra to construct a series of different individual PLSR models, and then aggregates those individual models.  The number of latent vectors used within each PLSR model is selected my minimizing a 10-fold average test set RMSE for increasing numbers of vectors (similar, though slightly different, than the A0 mode of the autoPLS package in R).

Example input and settings files are provided.  Input data is expected as a CSV with required columns for:
CSP_CODE [a unique identifier for each crown],
Band_# [all reflectance data for n bands included]
CalVal [indication of the calibration and validation subsets]
Chemistry values


If you use this code in a scientific publication, the authors would appreciate citations to the following publication:
R.E. Martin, K.D. Chadwick, P.G. Brodrick, L. Carranza-Jimenez, N.R. Vaughn, and G.P. Asner. An approach for foliar trait retrieval from airborne imaging spectroscopy of tropical forests. Remote Sensing (In Review), 2017.
