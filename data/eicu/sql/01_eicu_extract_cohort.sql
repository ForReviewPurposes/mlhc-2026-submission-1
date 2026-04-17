CREATE OR REPLACE TABLE `eicu-ext.eicu_ext_data.eicu_cohort_all` AS
WITH ranked AS (
  SELECT
    p.patientunitstayid,
    p.patienthealthsystemstayid,
    p.hospitalid,
    CASE
      WHEN p.age = '> 89' THEN 90
      WHEN p.age = '' THEN NULL
      ELSE CAST(p.age AS NUMERIC)
    END AS anchor_age,
    p.gender,
    p.hospitaladmitoffset,
    p.unitdischargeoffset,
    p.hospitaldischargeoffset,
    p.unitdischargestatus,
    p.hospitaldischargestatus,
    CASE WHEN p.hospitaldischargestatus = 'Expired' THEN 1 ELSE 0 END AS mortality,
    ROW_NUMBER() OVER (
      PARTITION BY p.patienthealthsystemstayid
      ORDER BY p.unitvisitnumber
    ) AS rn
  FROM `physionet-data.eicu_crd.patient` p
  WHERE p.hospitalid IN (73, 167, 264, 386, 328, 383)
    AND p.hospitaldischargestatus != ''
)
SELECT
  patientunitstayid,
  patienthealthsystemstayid,
  hospitalid,
  anchor_age,
  gender,
  hospitaladmitoffset,
  unitdischargeoffset,
  hospitaldischargeoffset,
  unitdischargestatus,
  hospitaldischargestatus,
  mortality,
  1440 AS landmark_offset,

  CASE
    WHEN hospitaldischargestatus = 'Expired'
         AND hospitaldischargeoffset >= 1440
    THEN 1
    ELSE 0
  END AS event_after_24h,

  ROUND((hospitaldischargeoffset - 1440) / 60.0, 0) AS duration_hours_from_24h,

  CASE
    WHEN hospitaldischargestatus = 'Expired'
         AND hospitaldischargeoffset >= 1440
    THEN ROUND((hospitaldischargeoffset - 1440) / 60.0, 0)
    ELSE NULL
  END AS time_to_death_from_24h_hours

FROM ranked
WHERE rn = 1
  AND anchor_age IS NOT NULL
  AND anchor_age >= 18
  AND hospitaldischargeoffset IS NOT NULL
  AND hospitaldischargeoffset >= 1440;