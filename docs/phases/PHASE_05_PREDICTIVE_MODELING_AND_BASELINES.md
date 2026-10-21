# Phase 5: Predictive Modeling & Baseline Evaluation

**Consolidated Steps**: Step 10 (Baseline Modeling), Step 10.1 (Independent Metrics Audit)  
**Primary Deliverables**: `reports/13_retrained_oof_predictions.csv`, `reports/10_baseline_metrics.csv`, `reports/10_family_level_metrics.csv`  

---

## 1. Scientific Objective
Evaluate whether protein representations (physicochemical properties vs. ESM-2 transformer embeddings) contain discriminative signal distinguishing pre-replicative regulatory enzymes ($\text{PRE} = \text{IE} + \text{E}$) from post-replicative structural virion proteins ($\text{POST} = \text{L1} + \text{L2}$) across unseen viral gene families.

## 2. Models Evaluated
1. **Majority Baseline**: `DummyClassifier(strategy='most_frequent')`
2. **Physicochemical Logistic Regression**: `LogisticRegression(C=1.0, penalty='l2', solver='lbfgs', max_iter=1000)` on 8 physical features.
3. **Physicochemical Random Forest**: `RandomForestClassifier(n_estimators=100, max_depth=5, min_samples_leaf=3)` on 8 physical features.
4. **Ridge-Regularized ESM-2**: `LogisticRegression(C=1.0, penalty='l2', solver='lbfgs', max_iter=1000)` on 480 ESM-2 dimensions.
5. **Combined Model**: `LogisticRegression(C=1.0, penalty='l2', solver='lbfgs', max_iter=1000)` on 488 concatenated features.

## 3. Evaluation Protocol
- **Generalization Unit**: Gene Family ($N=63$). Family out-of-fold probability is the mean sequence probability: $\bar{P}_{\text{PRE}} \ge 0.5 \implies \text{PRE}$.
- Preprocessing: `StandardScaler` fit strictly inside training folds.

## 4. Re-Verified Benchmark Results
| Model | Family Accuracy | Family BalAcc | Family Macro-F1 | Sequence Macro-F1 | Correct Families |
| :--- | :---: | :---: | :---: | :---: | :---: |
| Majority Baseline | 47.62% | 53.21% | 0.4741 | 0.2651 | 30 / 63 |
| Physicochemical RF | 58.73% | 51.79% | 0.5102 | 0.5399 | 37 / 63 |
| Physicochemical LogReg | 60.32% | 56.74% | 0.5679 | 0.5747 | 38 / 63 |
| Combined (Phys + ESM-2) | 77.78% | 76.96% | 0.7644 | 0.8966 | 49 / 63 |
| **Ridge ESM-2** | **80.95%** | **79.46%** | **0.7946** | **0.9112** | **51 / 63** |

## 5. Validation & Test Outcomes
- Validated under `tests/test_10_baseline_evaluation.py` and `test_10_1_independent_verification.py` (20 passing tests).
- Mathematical equality verified: For $N=63$ with 17 true PRE and 34 true POST correct, Balanced Accuracy $\equiv$ Macro-F1 $\equiv 0.7946$.
