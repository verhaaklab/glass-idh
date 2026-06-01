/*
GLASS-I: dNdScv Input Hypermutant Cohort
Author: C.M.S. Tesileanu
Date: 2026-05-27
*/
SELECT
  tp.case_barcode,
  (CASE WHEN pgeno.chrom = 23 THEN 'X' ELSE pgeno.chrom::varchar(2) END) AS chrom,
  lower(pgeno.pos) AS pos,
  pgeno.ref,
  pgeno.alt AS mut,
  st.idh_codel_subtype AS subtype,
    (CASE WHEN mutect2_call_a AND mutect2_call_b     THEN 'S'
          WHEN mutect2_call_a AND NOT mutect2_call_b THEN 'P'
          WHEN mutect2_call_b AND NOT mutect2_call_a THEN 'R' END) AS fraction
FROM variants.pgeno
INNER JOIN analysis.gold_set tp ON tp.tumor_pair_barcode = pgeno.tumor_pair_barcode
LEFT JOIN clinical.subtypes st ON st.case_barcode = pgeno.case_barcode
INNER JOIN analysis.tumor_clinical_comparison_2024 tcc ON tcc.tumor_pair_barcode = tp.tumor_pair_barcode
WHERE
    (mutect2_call_a OR mutect2_call_b) AND sbs11_r_hypermutation IS TRUE
    
-- END --