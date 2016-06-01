# weaklens_pipeline
Weak lensing pipeline for HSC

(1) The goal is to measure the HSC weak lensing signal around a set of different
lens catalogs, using a set of shape and photometric redshift catalogs from  HSC.

(2) Let us break this into steps: preprocess the lens and shape catalogs. In the
end as an input to weaklens.py, we need

    (a) For lens catalog:
        (i) ra, dec, z

    (b) For source catalog:
        (i) ra, dec, e1, e2, zmean, zmode, zmedian, zmc, z68up, z68down, z95up, z95down, and possibly p(z)
