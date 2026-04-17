CREATE OR REPLACE TABLE `mimic-iv-ext-icumort-1.mimiciv_ext_icumort_1.mimiciv_ext_icu_scores_2`
AS
WITH 
  scores AS (
    SELECT c.stay_id, a.apsiii, l.lods, o.oasis, sa.sapsii, si.sirs, so.sofa_24hours
    FROM `mimic-iv-ext-icumort-1.mimiciv_ext_icumort_1.mimiciv_ext_icu_cohort_1_timetodeath` c
    LEFT JOIN `physionet-data.mimiciv_3_1_derived.apsiii` a
    ON c.stay_id = a.stay_id
    LEFT JOIN `physionet-data.mimiciv_3_1_derived.lods` l
    ON c.stay_id = l.stay_id
    LEFT JOIN `physionet-data.mimiciv_3_1_derived.oasis` o
    ON c.stay_id = o.stay_id
    LEFT JOIN `physionet-data.mimiciv_3_1_derived.sapsii` sa
    ON c.stay_id = sa.stay_id
    LEFT JOIN `physionet-data.mimiciv_3_1_derived.sirs` si
    ON c.stay_id = si.stay_id
    LEFT JOIN `physionet-data.mimiciv_3_1_derived.sofa` so
    ON c.stay_id = so.stay_id
    WHERE so.hr=24
  )
SELECT * FROM scores
