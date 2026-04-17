CREATE OR REPLACE TABLE `mimic-iv-ext-icumort-1.mimiciv_ext_icumort_1.mimiciv_ext_icu_cohort_1_timetodeath` AS
WITH base AS (
  SELECT
    icu.subject_id,
    icu.hadm_id,
    icu.stay_id,
    icu.intime,
    icu.outtime,
    adm.dischtime,
    adm.deathtime,
    adm.hospital_expire_flag AS mortality,
    pat.gender,
    pat.anchor_age,
    adm.insurance,
    ROW_NUMBER() OVER (
      PARTITION BY icu.subject_id
      ORDER BY icu.intime
    ) AS rn
  FROM `physionet-data.mimiciv_3_1_icu.icustays` icu
  JOIN `physionet-data.mimiciv_3_1_hosp.admissions` adm
    ON icu.hadm_id = adm.hadm_id
  JOIN `physionet-data.mimiciv_3_1_hosp.patients` pat
    ON icu.subject_id = pat.subject_id
),
cohort AS (
  SELECT
    subject_id,
    hadm_id,
    stay_id,
    intime,
    outtime,
    dischtime,
    deathtime,
    mortality,
    gender,
    anchor_age,
    insurance,
    TIMESTAMP_ADD(intime, INTERVAL 24 HOUR) AS landmark_time
  FROM base
  WHERE rn = 1
    AND anchor_age >= 18
    AND dischtime IS NOT NULL
    AND dischtime >= TIMESTAMP_ADD(intime, INTERVAL 24 HOUR)
)
SELECT
  subject_id,
  hadm_id,
  stay_id,
  intime,
  outtime,
  dischtime,
  deathtime,
  mortality,
  gender,
  anchor_age,
  insurance,
  landmark_time,

  -- descriptive only
  CASE
    WHEN deathtime IS NOT NULL THEN TIMESTAMP_DIFF(deathtime, intime, HOUR)
    ELSE NULL
  END AS time_to_death_from_icu_hours,

  -- event after the 24h landmark
  CASE
    WHEN deathtime IS NOT NULL AND deathtime >= landmark_time THEN 1
    ELSE 0
  END AS event_after_24h,

  -- duration for survival analysis starting at the 24h landmark
  CASE
    WHEN deathtime IS NOT NULL AND deathtime >= landmark_time
      THEN TIMESTAMP_DIFF(deathtime, landmark_time, HOUR)
    ELSE TIMESTAMP_DIFF(dischtime, landmark_time, HOUR)
  END AS duration_hours_from_24h,

  -- optional explicit components
  CASE
    WHEN deathtime IS NOT NULL AND deathtime >= landmark_time
      THEN TIMESTAMP_DIFF(deathtime, landmark_time, HOUR)
    ELSE NULL
  END AS time_to_death_from_24h_hours,

  CASE
    WHEN dischtime >= landmark_time
      THEN TIMESTAMP_DIFF(dischtime, landmark_time, HOUR)
    ELSE NULL
  END AS time_to_hospital_discharge_from_24h_hours

FROM cohort
WHERE
  -- exclude deaths before the 24h landmark
  (deathtime IS NULL OR deathtime >= landmark_time);