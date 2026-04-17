CREATE OR REPLACE TABLE `eicu-ext.eicu_ext_data.eicu_vals_raw` AS

WITH cohort AS (
  SELECT * FROM `eicu-ext.eicu_ext_data.eicu_cohort_all`
),

vitals AS (
  SELECT
    c.patientunitstayid,
    c.hospitalid,
    v.chartoffset AS offset,
    v.heartrate AS heart_rate,
    v.respiratoryrate AS resp_rate,
    v.spo2,
    v.nibp_systolic AS sbp,
    v.nibp_diastolic AS dbp,
    v.nibp_mean AS mbp,
    v.temperature
  FROM cohort c
  JOIN `physionet-data.eicu_crd_derived.pivoted_vital` v
   ON c.patientunitstayid = v.patientunitstayid
  WHERE v.chartoffset BETWEEN 0 AND 1440

),

labs AS (
  SELECT
    c.patientunitstayid,
    c.hospitalid,
    l.chartoffset AS offset,
    l.creatinine,
    l.sodium,
    l.potassium,
    l.glucose,
    l.lactate,
    l.bicarbonate,
    l.wbc,
    l.platelets AS platelet,
    l.hemoglobin
  FROM cohort c
  JOIN `physionet-data.eicu_crd_derived.pivoted_lab` l
    ON c.patientunitstayid = l.patientunitstayid
  WHERE l.chartoffset BETWEEN 0 AND 1440
),

gcs_all AS (
  SELECT
    c.patientunitstayid,
    c.hospitalid,
    g.chartoffset AS offset,
    g.gcs,
    g.gcsmotor AS gcs_motor,
    g.gcsverbal AS gcs_verbal,
    g.gcseyes AS gcs_eyes
  FROM cohort c
  JOIN `physionet-data.eicu_crd_derived.pivoted_gcs` g
    ON c.patientunitstayid = g.patientunitstayid
  WHERE g.chartoffset BETWEEN 0 AND 1440

),

long AS (

SELECT patientunitstayid, hospitalid, offset, 'heart_rate' AS signal, heart_rate AS value FROM vitals WHERE heart_rate IS NOT NULL
UNION ALL

SELECT patientunitstayid, hospitalid, offset, 'mbp' AS signal, mbp AS value FROM vitals WHERE mbp IS NOT NULL
UNION ALL

SELECT patientunitstayid, hospitalid, offset, 'sbp' AS signal, sbp AS value FROM vitals WHERE sbp IS NOT NULL
UNION ALL

SELECT patientunitstayid, hospitalid, offset, 'dbp' AS signal, dbp AS value FROM vitals WHERE dbp IS NOT NULL
UNION ALL

SELECT patientunitstayid, hospitalid, offset, 'resp_rate' AS signal, resp_rate AS value FROM vitals WHERE resp_rate IS NOT NULL
UNION ALL

SELECT patientunitstayid, hospitalid, offset, 'spo2' AS signal, spo2 AS value FROM vitals WHERE spo2 IS NOT NULL
UNION ALL

SELECT patientunitstayid, hospitalid, offset, 'temperature' AS signal, temperature AS value FROM vitals WHERE temperature IS NOT NULL
UNION ALL

SELECT patientunitstayid, hospitalid, offset, 'creatinine' AS signal, creatinine AS value FROM labs WHERE creatinine IS NOT NULL
UNION ALL

SELECT patientunitstayid, hospitalid, offset, 'sodium' AS signal, sodium AS value FROM labs WHERE sodium IS NOT NULL
UNION ALL

SELECT patientunitstayid, hospitalid, offset, 'potassium' AS signal, potassium AS value FROM labs WHERE potassium IS NOT NULL
UNION ALL

SELECT patientunitstayid, hospitalid, offset, 'bicarbonate' AS signal, bicarbonate AS value FROM labs WHERE bicarbonate IS NOT NULL
UNION ALL

SELECT patientunitstayid, hospitalid, offset, 'hemoglobin' AS signal, hemoglobin AS value FROM labs WHERE hemoglobin IS NOT NULL
UNION ALL

SELECT patientunitstayid, hospitalid, offset, 'glucose' AS signal, glucose AS value FROM labs WHERE glucose IS NOT NULL
UNION ALL

SELECT patientunitstayid, hospitalid, offset, 'wbc' AS signal, wbc AS value FROM labs WHERE wbc IS NOT NULL
UNION ALL

SELECT patientunitstayid, hospitalid, offset, 'platelet' AS signal, platelet AS value FROM labs WHERE platelet IS NOT NULL
UNION ALL

SELECT patientunitstayid, hospitalid, offset, 'lactate' AS signal, lactate AS value FROM labs WHERE lactate IS NOT NULL
UNION ALL

SELECT patientunitstayid, hospitalid, offset, 'gcs' AS signal, gcs AS value FROM gcs_all WHERE gcs IS NOT NULL
UNION ALL

SELECT patientunitstayid, hospitalid, offset, 'gcs_motor' AS signal, gcs_motor AS value FROM gcs_all WHERE gcs_motor IS NOT NULL
UNION ALL

SELECT patientunitstayid, hospitalid, offset, 'gcs_verbal' AS signal, gcs_verbal AS value FROM gcs_all WHERE gcs_verbal IS NOT NULL
UNION ALL

SELECT patientunitstayid, hospitalid, offset, 'gcs_eyes' AS signal, gcs_eyes AS value FROM gcs_all WHERE gcs_eyes IS NOT NULL

)

SELECT *
FROM long;