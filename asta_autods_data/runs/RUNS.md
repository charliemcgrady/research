# AutoDiscovery Runs — Driver-Negative LUAD: CIN → immune escape (TIGIT/LAG3)

Hosted AutoDiscovery runs testing the theorizer's hypothesis. All 5 submitted
2026-07-23, **20 experiments each** (100 credits total), on the shared analysis-ready
table `derived/master_sample_table.csv` (501 TCGA-LUAD samples; 118 driver-negative /
383 driver-positive, RPA wild-type definition).

| Run | Gap | Hypothesis / intent | Run ID |
|-----|-----|---------------------|--------|
| 1 association | G1 | CIN burden ↔ TIGIT/LAG3 within driver-negative LUAD | `c9508770-af22-4441-aff0-ac028c2e831d` |
| 2 specificity | G2 | CIN×driver-status interaction — is the effect specific to driver-negative? | `34d6418f-0341-45b6-8d5f-93284984d592` |
| 3 selectivity | G3 | CIN elevates TIGIT/LAG3 selectively vs PD-1/PD-L1/TIM-3/CTLA-4? | `14653a15-4d45-4450-8a90-b62a2c61d61f` |
| 4 mechanism | G4 | Is CIN→TIGIT/LAG3 mediated by IFN/cGAS-STING? | `8ceb95fa-2cf1-413c-806c-aaef41111044` |
| 5 CN events | G5 | Which specific CN events (9p21/B2M/TP53/RB1/…) predict checkpoint-high? | `d9363c9e-b5b9-4b21-afc1-35045f79ba44` |

Config: runs 1–4 `exploration_weight=2.0, surprisal_width=0.4, evidence_weight=2.5`;
run 5 (exploratory) `4.0 / 0.3 / 2.0`. All `mcts_selection=ucb1_recursive`.
Every run's intent instructs adjusting for **TMB** (confounder).

## Monitoring
```bash
export PATH="$HOME/.local/bin:$PATH"
asta autodiscovery status <runid>
asta autodiscovery experiments <runid> --format text
asta autodiscovery run <runid> --format text
```
