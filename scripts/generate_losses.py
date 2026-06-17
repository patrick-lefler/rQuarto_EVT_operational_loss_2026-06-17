"""
generate_losses.py
------------------
Synthetic operational loss register for NexaCore Financial Technologies.
Produces two CSV outputs:
  - data/loss_register.csv         (~500 loss events, 2020-2024)
  - data/threshold_diagnostics.csv (mean excess plot values across candidate thresholds)

Architecture
------------
Loss severity is a two-component mixture:
  Bulk component  (~85% of events): log-normal, calibrated to everyday operational losses
  Tail component  (~15% of events): Pareto, calibrated to produce GPD shape xi ~ 0.3-0.4

This mixture ensures the log-normal fit to the full register visibly underfits the right
tail — the central visual argument of the document.

Basel II event type frequency allocation:
  Execution / Delivery / Process Management : 35%
  External Fraud                             : 20%
  Clients / Products / Business Practices   : 15%
  Business Disruption / System Failures     : 12%
  Internal Fraud                             :  8%
  Employment Practices                       :  6%
  Damage to Physical Assets                  :  4%

All randomness seeded at 42 for reproducibility.

Author : Patrick Lefler
Project: The Tail You Cannot Fit — EVT for Operational Loss Capital Estimation
"""

import os
import numpy as np
import pandas as pd
from scipy import stats
from datetime import date, timedelta

# ---------------------------------------------------------------------------
# 0. Seed and output directory
# ---------------------------------------------------------------------------
RNG = np.random.default_rng(42)
os.makedirs("data", exist_ok=True)

# ---------------------------------------------------------------------------
# 1. Configuration
# ---------------------------------------------------------------------------
N_EVENTS       = 500
DATE_START     = date(2020, 1, 2)
DATE_END       = date(2024, 12, 31)

# Mixture proportions
BULK_FRAC      = 0.85
TAIL_FRAC      = 0.15

# Bulk component: log-normal
# ln(X) ~ N(mu, sigma) => median = exp(mu), mean = exp(mu + sigma^2/2)
# Target: median gross loss ~ EUR 12,000; moderate right skew
BULK_LOGNORM_MU    = np.log(12_000)   # log-scale mean  => median EUR 12,000
BULK_LOGNORM_SIGMA = 1.10             # log-scale SD    => mean  ~ EUR 27,000

# Tail component: Pareto (scipy parameterisation: Pareto(b) where b = 1/xi)
# Target xi ~ 0.35 => b = 1/0.35 ~ 2.857
# Scale (x_min) set so tail losses start above EUR 80,000
TAIL_PARETO_B      = 2.857            # shape (1/xi)
TAIL_PARETO_SCALE  = 80_000          # minimum tail loss (EUR)

# Recovery rate: bulk losses recover 5-25%, tail losses recover 0-10%
BULK_RECOVERY_LO,  BULK_RECOVERY_HI  = 0.05, 0.25
TAIL_RECOVERY_LO,  TAIL_RECOVERY_HI  = 0.00, 0.10

# Basel II event types and their frequency weights
BASEL_TYPES = [
    "Execution / Delivery / Process Management",
    "External Fraud",
    "Clients / Products / Business Practices",
    "Business Disruption / System Failures",
    "Internal Fraud",
    "Employment Practices",
    "Damage to Physical Assets",
]
BASEL_WEIGHTS = np.array([0.35, 0.20, 0.15, 0.12, 0.08, 0.06, 0.04])

# NexaCore business lines (plausible for an EU fintech)
BUSINESS_LINES = [
    "Payments & Settlement",
    "Retail Banking",
    "Corporate Banking",
    "Asset Management",
    "Technology & Operations",
    "Compliance & Legal",
]
BL_WEIGHTS = np.array([0.30, 0.25, 0.18, 0.12, 0.10, 0.05])

# Tail event type over-representation: tail losses skew toward these categories
TAIL_BASEL_OVERWEIGHTS = {
    "External Fraud":                            0.35,
    "Clients / Products / Business Practices":   0.30,
    "Business Disruption / System Failures":     0.20,
    "Execution / Delivery / Process Management": 0.10,
    "Internal Fraud":                            0.05,
    "Employment Practices":                      0.00,
    "Damage to Physical Assets":                 0.00,
}
TAIL_BASEL_WEIGHTS = np.array([TAIL_BASEL_OVERWEIGHTS[t] for t in BASEL_TYPES])

# ---------------------------------------------------------------------------
# 2. Generate event dates (uniformly distributed over range, then sorted)
# ---------------------------------------------------------------------------
total_days = (DATE_END - DATE_START).days
day_offsets = np.sort(RNG.integers(0, total_days, size=N_EVENTS))
event_dates = [DATE_START + timedelta(days=int(d)) for d in day_offsets]

# ---------------------------------------------------------------------------
# 3. Assign component membership
# ---------------------------------------------------------------------------
n_bulk = int(np.round(N_EVENTS * BULK_FRAC))
n_tail = N_EVENTS - n_bulk

component = np.array(["bulk"] * n_bulk + ["tail"] * n_tail)
RNG.shuffle(component)

# ---------------------------------------------------------------------------
# 4. Sample gross losses
# ---------------------------------------------------------------------------
gross_losses = np.empty(N_EVENTS)

bulk_idx = np.where(component == "bulk")[0]
tail_idx = np.where(component == "tail")[0]

# Bulk: log-normal
bulk_raw = RNG.lognormal(mean=BULK_LOGNORM_MU, sigma=BULK_LOGNORM_SIGMA, size=len(bulk_idx))
# Floor at EUR 500 (minimum recordable loss)
bulk_raw = np.maximum(bulk_raw, 500.0)
gross_losses[bulk_idx] = np.round(bulk_raw, 2)

# Tail: Pareto (scipy: pareto.rvs(b, scale=scale) => X ~ scale * (1 + Pareto(b))
# scipy.stats.pareto(b, loc=0, scale=s) gives X >= s with P(X>x) = (s/x)^b
tail_raw = stats.pareto.rvs(
    b=TAIL_PARETO_B,
    scale=TAIL_PARETO_SCALE,
    size=len(tail_idx),
    random_state=np.random.RandomState(seed=99)   # deterministic secondary seed
)
gross_losses[tail_idx] = np.round(tail_raw, 2)

# ---------------------------------------------------------------------------
# 5. Sample recovery amounts
# ---------------------------------------------------------------------------
recovery_rates = np.empty(N_EVENTS)

bulk_rec_rates = RNG.uniform(BULK_RECOVERY_LO, BULK_RECOVERY_HI, size=len(bulk_idx))
tail_rec_rates = RNG.uniform(TAIL_RECOVERY_LO, TAIL_RECOVERY_HI, size=len(tail_idx))

recovery_rates[bulk_idx] = bulk_rec_rates
recovery_rates[tail_idx] = tail_rec_rates

recovery_amounts = np.round(gross_losses * recovery_rates, 2)
net_losses       = np.round(gross_losses - recovery_amounts, 2)

# ---------------------------------------------------------------------------
# 6. Assign Basel event types
# ---------------------------------------------------------------------------
basel_bulk = RNG.choice(BASEL_TYPES, size=len(bulk_idx), p=BASEL_WEIGHTS / BASEL_WEIGHTS.sum())
basel_tail = RNG.choice(BASEL_TYPES, size=len(tail_idx), p=TAIL_BASEL_WEIGHTS / TAIL_BASEL_WEIGHTS.sum())

basel_event_type = np.empty(N_EVENTS, dtype=object)
basel_event_type[bulk_idx] = basel_bulk
basel_event_type[tail_idx] = basel_tail

# ---------------------------------------------------------------------------
# 7. Assign business lines
# ---------------------------------------------------------------------------
business_line = RNG.choice(BUSINESS_LINES, size=N_EVENTS, p=BL_WEIGHTS / BL_WEIGHTS.sum())

# ---------------------------------------------------------------------------
# 8. Assemble loss register
# ---------------------------------------------------------------------------
loss_register = pd.DataFrame({
    "event_date":      [d.isoformat() for d in event_dates],
    "business_line":   business_line,
    "basel_event_type": basel_event_type,
    "gross_loss_eur":  gross_losses,
    "recovery_eur":    recovery_amounts,
    "net_loss_eur":    net_losses,
    "component":       component,   # kept for validation; not used in document
})

# ---------------------------------------------------------------------------
# 9. Validation checks
# ---------------------------------------------------------------------------
print("=" * 60)
print("LOSS REGISTER VALIDATION")
print("=" * 60)
print(f"Total events          : {len(loss_register)}")
print(f"Bulk events           : {(loss_register.component == 'bulk').sum()}")
print(f"Tail events           : {(loss_register.component == 'tail').sum()}")
print(f"Date range            : {loss_register.event_date.min()} to {loss_register.event_date.max()}")
print()
print(f"Gross loss (EUR)")
print(f"  Min                 : {loss_register.gross_loss_eur.min():>12,.0f}")
print(f"  Median              : {loss_register.gross_loss_eur.median():>12,.0f}")
print(f"  Mean                : {loss_register.gross_loss_eur.mean():>12,.0f}")
print(f"  95th pctile         : {loss_register.gross_loss_eur.quantile(0.95):>12,.0f}")
print(f"  99th pctile         : {loss_register.gross_loss_eur.quantile(0.99):>12,.0f}")
print(f"  Max                 : {loss_register.gross_loss_eur.max():>12,.0f}")
print()
print(f"Net loss (EUR)")
print(f"  Min                 : {loss_register.net_loss_eur.min():>12,.0f}")
print(f"  Median              : {loss_register.net_loss_eur.median():>12,.0f}")
print(f"  Mean                : {loss_register.net_loss_eur.mean():>12,.0f}")
print(f"  95th pctile         : {loss_register.net_loss_eur.quantile(0.95):>12,.0f}")
print(f"  99th pctile         : {loss_register.net_loss_eur.quantile(0.99):>12,.0f}")
print(f"  Max                 : {loss_register.net_loss_eur.max():>12,.0f}")
print()
print("Basel event type distribution:")
print(loss_register.groupby("basel_event_type")["gross_loss_eur"]
      .agg(count="count", median="median", mean="mean")
      .round(0).to_string())

# ---------------------------------------------------------------------------
# 10. Threshold diagnostics (mean excess plot data)
# ---------------------------------------------------------------------------
# For threshold u, mean excess = E[X - u | X > u]
# Compute over a grid of candidate thresholds from 10th to 97th percentile

net = loss_register.net_loss_eur.values
thresholds = np.percentile(net, np.arange(10, 97, 1))
thresholds = np.unique(np.round(thresholds, 0))

mep_rows = []
for u in thresholds:
    exceedances = net[net > u] - u
    n_exceed = len(exceedances)
    if n_exceed < 5:
        break
    mean_excess = exceedances.mean()
    # 95% CI via CLT approximation
    se = exceedances.std(ddof=1) / np.sqrt(n_exceed)
    mep_rows.append({
        "threshold_eur":    float(u),
        "n_exceedances":    int(n_exceed),
        "mean_excess_eur":  float(np.round(mean_excess, 2)),
        "ci_lower_eur":     float(np.round(mean_excess - 1.96 * se, 2)),
        "ci_upper_eur":     float(np.round(mean_excess + 1.96 * se, 2)),
    })

threshold_diagnostics = pd.DataFrame(mep_rows)

print()
print("=" * 60)
print("THRESHOLD DIAGNOSTICS")
print("=" * 60)
print(f"Threshold grid        : {len(threshold_diagnostics)} candidate values")
print(f"Threshold range       : EUR {threshold_diagnostics.threshold_eur.min():,.0f} "
      f"to EUR {threshold_diagnostics.threshold_eur.max():,.0f}")
# Flag approximate linearity zone (where mean excess is still rising, indicating GPD validity)
# Heuristic: look for the stretch where mean excess gradient is positive and relatively stable
me_vals = threshold_diagnostics.mean_excess_eur.values
me_thresh = threshold_diagnostics.threshold_eur.values
# Exceedance counts at key round thresholds for document reference
for u_check in [50_000, 75_000, 100_000, 125_000, 150_000]:
    row = threshold_diagnostics[threshold_diagnostics.threshold_eur <= u_check]
    if len(row):
        r = row.iloc[-1]
        print(f"  u = EUR {u_check:>7,.0f}  =>  n_exceedances = {r.n_exceedances:>4}  "
              f"mean_excess = EUR {r.mean_excess_eur:>10,.0f}")

# ---------------------------------------------------------------------------
# 11. Write outputs
# ---------------------------------------------------------------------------
loss_register.to_csv("data/loss_register.csv", index=False)
threshold_diagnostics.to_csv("data/threshold_diagnostics.csv", index=False)

print()
print("=" * 60)
print("OUTPUT FILES WRITTEN")
print("=" * 60)
print("  data/loss_register.csv         rows:", len(loss_register))
print("  data/threshold_diagnostics.csv rows:", len(threshold_diagnostics))
print()
print("Next step: run compute_evt.R to fit GPD and produce gpd_fit.csv")
print("           and capital_comparison.csv")
