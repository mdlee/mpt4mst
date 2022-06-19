The `.mat` file is a MATLAB data file, with a single structured variable `d` and the following fields, using a long format for the trial-by-trial behavior:

```
                  trialON: [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 … ]
                lureBinON: [3 2 5 4 5 2 5 1 5 1 5 2 2 3 1 1 2 4 2 1 4 5 2 4 1 2 … ]
                   lureON: [1 0 0 1 0 0 0 0 0 1 0 1 0 0 0 0 1 0 0 1 0 0 0 0 0 0 … ]
          participantIDON: [47878 47878 47878 47878 47878 47878 47878 47878 … ]
                  truthON: [2 2 1 2 1 1 2 1 2 2 2 2 2 2 1 1 2 1 1 2 2 1 1 1 2 1 … ]
                  studyON: [0 0 42 0 23 89 0 61 0 0 0 0 0 0 38 43 0 29 65 0 0 … ]
                    gapON: [129 130 89 132 110 45 135 75 137 138 139 140 141 … ]
                 gapIdxON: [35 36 19 38 27 5 41 14 43 44 45 46 47 48 24 22 50 … ]
               decisionON: [NaN 2 1 2 1 1 2 1 2 1 2 1 2 2 1 2 1 1 1 1 2 1 2 1 2 … ]
                correctON: [0 1 1 1 1 1 1 1 1 0 1 0 1 1 1 0 0 1 1 0 1 1 0 1 1 1 … ]
            participantON: [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 … ]
               stimulusON: [246 334 234 284 217 240 311 263 362 312 224 227 243 … ]
     participantCorrectON: [150 165 145 94 153 150 125 149 145 144 129 137 147 … ]
                nTrialsON: [192 192 192 192 192 192 192 192 192 192 192 192 192 … ]
                 trialOSN: [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 … ]
               lureBinOSN: [3 2 5 5 1 5 3 3 5 4 1 5 1 5 4 4 2 3 3 3 3 5 1 3 4 5 … ]
                  lureOSN: [1 1 0 0 1 0 1 0 0 0 0 0 0 1 0 0 0 0 0 0 1 0 1 1 1 0 … ]
         participantIDOSN: [47878 47878 47878 47878 47878 47878 47878 47878 … ]
                 truthOSN: [3 3 1 2 3 1 3 2 1 2 2 1 1 3 2 1 1 2 2 2 3 2 3 3 3 1 … ]
                 studyOSN: [118 37 126 0 62 58 13 0 119 0 0 47 120 59 0 80 87 0 … ]
                   gapOSN: [11 93 5 132 71 76 122 136 18 138 139 93 21 83 143 … ]
                gapIdxOSN: [2 41 1 73 26 28 64 76 3 78 79 41 4 33 83 21 17 86 … ]
              decisionOSN: [NaN 3 NaN 2 3 1 3 3 1 3 NaN 2 1 3 2 1 1 2 2 2 2 2 … ]
               correctOSN: [0 1 0 1 1 1 1 0 1 0 0 0 1 1 1 1 1 1 1 1 0 1 0 1 1 1 … ]
           participantOSN: [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 … ]
              stimulusOSN: [28 151 2 44 174 207 94 82 89 152 107 153 18 49 99 … ]
    participantCorrectOSN: [81 144 124 80 134 136 116 145 115 131 114 126 157 … ]
               nTrialsOSN: [192 192 192 192 192 192 192 192 192 192 192 192 192 … ]
     uniqueParticipantIDs: [47878 48757 54031 60382 61738 63091 63505 65224 … ]
            nParticipants: 21
                   nLures: 5
            studyTrialsON: 128
           studyTrialsOSN: 128
             uniqueGapsON: [5 18 21 38 45 47 53 57 58 64 66 69 70 75 76 82 83 … ]
            uniqueGapsOSN: [5 11 18 21 24 29 33 34 38 39 44 45 47 52 53 57 58 … ]
                  nGapsON: 209
                 nGapsOSN: 486
                      LDI: [0.1250 0.6406 0.0625 -0.0312 0.3906 0.2500 0.1094 … ]
                      REC: [0.4531 0.6562 0.9219 0.2969 0.7344 0.8750 0.7031 … ]
                     dpTF: [3.1961 2.9839 2.8849 -0.3297 2.6754 3.7158 1.9493 … ]
                     dpTL: [0.0760 0.5596 -0.2691 0.3894 0.3488 -0.1911 -0.5596 … ]
  ```

Core data. (+ indicates required input to JAGS model)
- `participantON`: the participant number who completed this trial
- `trialON`: the trial number
- `lureON`+: 0 = not a lure item, 1 = is a lure item
- `lureBinON`+: lure bin 1, 2, 3, 4, or 5.
- `decisionON`+: observed response, nan = missing, 1 = old, 2 = new
- `truthON`+: true nature of item, 1 = old, 2 = new. Note that the MATLAB code creates a variable `truth` that modifies `truthON` to be 3 for lure items.
- `correctON`: 1 = response was correct, 0 = incorrect

Extra stuff:
- `studyON`: the position in the study sequence that the currently test item was originally presented
- `gapON`: the number of trials between original study and current test for item
- 
