# Numerical Simulations for a Two-Public Good Agreement Game

This markdown is prepared to accompany the Working Paper "Global Governance on Multiple Public Goods". Only necessary information is repeated for completeness and understanding of the code. Please refer to the paper before reading this file. Glossary.md summarises the key terms used in this markdown.
 
We implement numerical simulations in Python, based on equilibrium solutions calculated in the mathematical software, Maple. 

The specific functions analysed is:

$\pi_i=aQ-\frac{\alpha}{2}Q^2+bX-\frac{\beta}{2}X^2+\gamma QX- \frac{m}{q_i}q_i^2 - \frac{g}{2}x_i^2$

The following steps are taken for the numerical simulations:

## Step 1. Transfer and Store the Mathematical Functions in Python

We simulate two versions of the game:

1. **Full Agreement (FA)** – where the coalition cooperates on both public goods.
2. **Q Agreement (QA)** – where cooperation is only on Q.

The mathematical results (equilibrium levels, stability, payoffs, etc.) are first solved in *Maple* and then converted to Python code by using the following command in Maple:

```
CodeGeneration[Python](function_name)
```

This gives us the following functions for each agreement type:

```
FA: faQ, faX, sigma_q, sigma_x, faQCOH, faXCOH, faPEP, faPIP, faSAD, faWCOH, faIS, faINI  
QA: qaQ, qaX, sigma_q, sigma_x, qaQCOH, qaXCOH, qaPEP, qaPIP, qaSAD, qaWCOH, qaIS, qaINI
```

> All functions are defined in the "Glossary.md". 

Which are copied from Maple then are stored in the files:
- `FA_functions.py` for Full Agreement  
- `QA_functions.py` for Q Agreement

---

### Parameters

To keep the code clean and consistent, we use a separate file called `parameters.py`. It groups all parameters into a single object:

```python
class Parameters:
    def __init__(self, a, b, m, g, alpha, beta, gamma, s):
        self.a = a
        self.b = b
        self.m = m
        self.g = g
        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma
        self.s = s
```

All functions use this structure so the same parameters are passed in the same way every time.

---

### Example Function

Each function follows the same format:

```python
def faQ(params: Parameters) -> float:
    a = params.a
    b = params.b
    ...
    return some_expression  # from Maple
```

This means:
- The function uses the values in the `Parameters` object.
- It returns a number (a float).
- The expression is the direct result from Maple.


## Step 2. Define Restrictions in Python

To ensure that our simulations are meaningful, we need to apply **restrictions** on the parameter values, as defined in Appendix X. These restrictions make sure that the solutions we analyze are **valid**, such as being **interior** (not at boundaries) (Cond...) and **unique** (only one equilibrium) (Cond...).

Because these restrictions are logically separate from the mathematical functions, we store them in a separate file: `restrictions_groups.py`.

---

### Types of Restrictions

Each game version (FA and QA) has **three key restrictions**:

- Two for ensuring interior solutions.
- One for ensuring uniqueness.

These restrictions are upper bounds on `gamma squared` and are labeled:

```
FA: FA_intq, FA_intx, FA_rest  
QA: QA_intq, QA_intx, QA_rest
```

Depending on the **sign of the cross-effects** (captured by `gamma`), we split the restrictions into two groups:

- **Positive cross-effects** (γ > 0):  
  - FA: `fa_pxr = [FA_rest]`  
  - QA: `qa_pxr = [QA_rest]`

- **Negative cross-effects** (γ < 0):  
  - FA: `fa_nxr = [FA_intq, FA_intx, FA_rest]`  
  - QA: `qa_nxr = [QA_intq, QA_intx, QA_rest]`

As the interior solutions always hold for positive cross effects, it can be removed from the restriction lists. Thus, the restrictions applied on gamma squared in the simulations, are cross effect specific restrictions.

---

### Sigma Groups

We also want to categorize parameter sets based on the **sign of the reaction functions** (`sigma_Q`, `sigma_X`). These signs depend on γ². 

Technical details as derived in the Maple file section "Sigma groups":
> Sigma is increasing in $s$, so it will only change signs once, going from negative to positive.
> 
> If we are interested in parameter sets in which sigma is negative for all $s$, then we impose $s=n$ such that this condition holds for all $s$. Doing so creates an upper bound on γ², as with the interior and uniqueness conditions.
> 
> If we are interesting in parameter sets in which sigma is positive for all $s$, we impose $s=1$. A lower bound on γ² is created, which, when combined with the uniqueness & interior conditions, will give us a relevant interval for γ².
> 
> As we have a sigma for both X and Q between-player reaction functions, we can have up to two lower bounds, and two (additional) upper bounds on γ².  

There are nine sigma groups:

- **G1_1**: Sigmas Q can change sign. Sigma X always positive.
- **G1_2**: Sigmas Q can change sign. Sigma X always negative.
- **G1_3**: Sigmas X can change sign. Sigma Q always positive.
- **G1_4**: Sigmas X can change sign. Sigma Q always negative.
- **G1_5(referred to as G1 in the code)**: Both sigmas can change sign (from negative to positive). 
- **G2**: Both sigmas are negative.
- **G3**: σ_Q is positive, σ_X is negative.
- **G4**: σ_Q is negative, σ_X is positive.
- **G5**: Both sigmas are positive.

Side notes:
> When sigma changes signs, we have both a lower and upper bound on γ².
> 
> In the Q agreement, Sigma_Q is not a function of $s$, thus G1_1, G1_2, and G1_5 are not relevant.  
---

### Sigma Restriction Functions

These functions define bounds on γ² which dictate the the **sign** of the reaction functions (`sigma_Q`, `sigma_X`). For example:

```python
# σ_Q always negative (FA and QA different)
def FA_sqn(a, b, g, d, alpha, beta, m, n):
    return alpha * (g + beta * n**2) / n**2

def QA_sqn(a, b, g, d, alpha, beta, m, n):
    return alpha * (beta * n + g) / n
```

```python
# σ_Q always positive (FA and QA the same)
def sqp(a, b, g, d, alpha, beta, m, n):
    return alpha * (beta * n + g) / n
```

- `sqn` gives the **upper bound** on γ² to ensure σ_Q remains **negative**.  
- `sqp` gives the **lower bound** on γ² to ensure σ_Q remains **positive**.

> The same structure applies to σ_X via `sxn` and `sxp`.
>  
> Functions with a `_c` suffix (e.g. `sqn_c`) handle **sign-change cases**—i.e. when a sigma changes sign as `s` increases. These create both a lower and upper bound on γ², defining a valid range rather than a single condition.

The sigma restriction functions are used to define the sigma groups defined above.
- `sqn` and/or `sxn` are listed in `G[1.1,1.2, etc]max`;
- `sqp` and/or `sxp` are listed in `G[1.1,1.2, etc]min`. 

---

### Practical Limits

To avoid computational issues:
- When sigma is **positive**, we keep `alpha` and `beta` **very small** (e.g. `0-2e-8`, which is defined as `max_sp`). 
- When sigma is **negative**, a full uniform range (defined below in normalisation) is allowed.

---

### Group Definitions

At the end of `restrictions_groups.py`, we group all the relevant restrictions together for the main file to distinguish different groups. Example:

```python
FA_group_definitions = {
    "G1": {
        "max_funcs": G1max,
        "min_funcs": G1min,
        "max_alpha": max_sp,
        "max_beta": max_sp
    },
    ...
}
```

This allows the simulation code below to loop through groups in a general way without hardcoding different bounds for different groups within the simulation - rather the group definition is "called" and the simulation is conducted for each group. 


## Step 3. Random Sampling of Parameters

We now generate random values for the model's parameters and check which combinations are valid. We predefine `n=200`.


#### Normalisation

We normalise the parameter range between 0 and 1 based on the observation that uniform scaling of parameter sets does not affect the qualitative outcomes of the model.

---

There are **two layers** of iteration for the sampling process. The first tests different forms of the benefit function. The second is the internal loop which focuses on sigma-group specific parameter sets. 


### External Loop: Testing Functional Forms

The **external loop** (defined by `rounds_config`) goes through four different rounds in order to test the different functional forms of the specific function:

- **Round I**: Both `alpha` and `beta` > 0 (strictly concave benefits for both public goods)  
- **Round II**: `alpha = 0`, `beta > 0` (linear benefits for Q)  
- **Round III**: `alpha > 0`, `beta = 0` (linear benefits for X)  
- **Round IV**: `alpha = 0`, `beta = 0` (linear benefits for both)  

> Why have separate rounds? This allows us to observe if the linear benefits fundamentally change any results. If 0 was simply included in the random sampling, we would not see directly if this plays a role.
>
> Note that if benefits are linear, the associated sigma will always be positive. 


### Internal Loop: Sampling and Evaluating Parameters

Within each round, we loop through each **sigma group** (`G1_1`, `G2`, etc.). Each group has its own restriction bracket on γ², and random parameters are sampled within these bounds. This ensures that all valid samples generated under a given group definition will by construction belong to that group.

> Note! A parameter set is not assigned to a group _after_ sampling; rather, the group definitions determine which parameter sets are admissible in the first place. This allows us to pre-sort parameter sets, and ensure that each group has an equal number of examples from which we can interpret. 

Within each group there are four steps:

**3.1.** **`(a, b)`** and **`(m, g)`** sampled from $\sim$ `Uniform (0, 1)`

**3.2.** **` (alpha, beta) `** sampled from $\sim$ `Uniform (0, bound)`  

**`bound` in round I**
   - As mentioned in the section on practical limits, when sigma is positive, it is necessary to implement an upper bound, `bound=max_sp` due to computational limitations.
   - If the sigmas are negative, `bound = 1`.  

**`bound` in rounds II-IV**
   - For linear benefits in Q: `alpha = 0`; For linear benefits in X: `beta = 0`
   - Non-linear benefits with negative sigma: `bound = 1`.
   - Non-linear benefits with positive sigma: `bound = max_sp`.

**3.3.** **Compute restriction bounds** for γ²:  
To ensure the parameter set will generate a unique, interior solution within the specified sigma group, γ² should be sampled from a specific range:
   -  `max_funcs` (different for each sigma group) combined with one of the lists below (applicable per agreement type and cross effect), define the upper bound of γ².
> `fa_nxr` for FA negative cross effects;  
> `fa_pxr` for FA positive cross effects;  
> `qa_nxr` for QA negative cross effects; or  
> `qa_pxr` for QA positive cross effects.
   - If at least one sigma is positive for some s, then there will also be a lower bound on γ², contained in `min_funcs`.

The relevant range for γ² is thus defined by a set of Upper bounds i.e. right hand side of inequality: `RHS = max_func + (e.g.)fa_nxr`, and lower bounds i.e. left hand side of inequality: `LHS = min_funcs`.

   - To find the most restrictive conditions, the randomly sampled parameters are substituted into the RHS and LHS, which will give a list of real numbers. The smallest number is retained as the smallest maximum of γ². The largest number is saved as the largest minimum of γ². 

> For a given parameter set, it is possible that the LHS (min of gamma sqrd) > RHS (max of gamma sqrd). Such parameter sets will be removed as they are invalid. To ensure an even set of parameters is sampled for each group, paramater sets continue to be sampled (Steps 3.1 - 3.3) until a pre-specified number of valid samples is reached. 

**3.4.** **Sample γ²** from within the valid range:  
   - **`(gamma_squared)`** sampled from $\sim$ `Uniform (LHS, RHS)`
   - Compute `gamma_pxr = +√γ²` or `gamma_nxr = -√γ²`, corresponding to **positive** and **negative** cross-effects.

We now have a full set of parameters, excluding 's'.


## Step 4. Compute objective function from 's=1..200'

In this step, we determine the sign of the internal stability function (PAPER reference) for each valid parameter set, for each `s=1..200`. 

**4.1.** **Loop over s** from 1 to 200:  
   - For each valid parameter set:  
     - Find the largest coalition size `s` such that the objective function (e.g., `faIS ≥ 0`) holds.  
     - Save results for each valid case.
> The largest stable coalition size is relevant as if s* is the largest stable coalition – any coalition above this is not internally stable. s* must also be externally stable. And we also know that for $\hat{s}<s^{*}$ which is also internally stable must be Pareto dominated according to the proof in Appendix X. 

## Step 5. Post processing: remaining functions and table compilation

For each valid parameter set, the following functions (mentioned in Step 1) are evaluated by iterating over 's=1..200'. 

```
FA: faQCOH, faXCOH, faPEP, faPIP, faSAD, faWCOH
QA: qaQCOH, qaXCOH, qaPEP, qaPIP, qaSAD, qaWCOH
```
The code records whether, for a given parameter set, the function is positive and/or negative for a given s, which can lead to the following outputs:
```
+... if positive for all s.
-... if negative for all s.
+/-... if positive for some s.
0... if zero for all s. 
```
The result for each function is then saved for the parameter set. 

The final function to evaluate is the improvement over no cooperation index `faINI, qaINI`. The largest stable coalition size found in step 4.1 is substituted into the relevant INI to determine the relative increase in welfare due to the stable coalition. The reported value is a percentage. 


**Table Compilation**  

Each parameter set generates a row in a table with the following information:
  - Round type (I–IV)
  - Group label (e.g. G2)
  - Sign of γ (positive or negative)
  - (in separate columns) the signs of `[]QCOH, []XCOH, []PEP, []PIP, []SAD, []WCOH` over $s$
  - The largest stable coalition size, $s*$
  - INI for the given $s*$.

A summary table is then generated to count the outcomes for a given round and group. That is, the no. of times QCOH>0 for all s, QCOH<0 for all s, etc occurred for a specific group within the round. 


  

  
---

## Appendix A: Technical Implementation Notes

The simulation script includes several implementation features to ensure robustness and reproducibility:

- **Invalid Range Filtering when there is a lower bound on γ²**  
  When the calculated bounds on γ² are not valid (e.g. LHS > RHS or non-finite values), those samples are discarded:
  ```python
  valid_mask = (RHS > LHS) & np.isfinite(LHS) & np.isfinite(RHS)
  ```
This function filters out many parameter sets. 

- **Separate Results for γ > 0 and γ < 0**  
  The code stores results for positive and negative γ values separately to allow for direct comparison:
  ```python
  valid_samples_fpos = [...]
  valid_samples_fneg = [...]
  ```

- **Objective Evaluation with Error Handling**  
  The objective functions (e.g. `faIS`) are evaluated inside a loop for each `s`, with exception handling for numerical issues:
  ```python
  try:
      val = faIS(params)
  except ZeroDivisionError:
      continue
  ```

- **Post-Check for Sigma Group Membership**  
  After storing valid samples, the code verifies that each sample still belongs to the correct sigma group:
  ```python
  run_sigma_checks(valid_samples_fpos, sigma_checks)
  ```

