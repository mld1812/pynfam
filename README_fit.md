# README for pynfam_fit_wrapper

## Overview
A python wrapper for fitting is provided in PyNFAM. 
The function `pynfam_fit_wrapper` constructs inputs for function `pynfam_mpi_calc`, launch the function and collect outputs from disk files. 
Four observables are supported by `pynfam_fit_wrapper`: Gamow-Teller resonance (GTR) energy, $\beta^-$ decay half life, spin-dipole resonance (SDR) energy and pairing gap. 
Only first two observables are needed for the fitting of time-odd Skyrme parameters and isoscalar pairing strength; 
the SDR energy is used to check our results, and the pairing gap is for the fitting of isovector pairing. 
In the current stage we are not going to apply advanced statistical tools on the fitting of isovector pairing.

`pynfam` is a Python library that provides inputs for Fortran executables, launch them and then collect outputs. 
There are three stages of one `pynfam` run: hfb, fam, beta. 
In the hfb stage `pynfam` runs Fortran program "HFBTHO" to obtain the ground state of nuclei. 
In the fam stage `pynfam` runs Fortran program "PNFAM" to calculate the linear response to various $\beta$ decay transitions. 
The hfb stage does not depend on the time-odd Skyrme parameters and isoscalar pairing strength, but the fam stage does depend on them. 
In the beta stage `pynfam` collects outputs of the fam stage, calculate the strength distributions and decay half lives (see "Summary of PNFAM calculation" for details), and write them to files. 

One `pynfam_mpi_calc` run can only deal with one type of observable, so `pynfam_fit_warpper` will sequentially launch several `pynfam_mpi_calc` runs if two or more observable types are wanted. 
For one MPI-parallelized `pynfam_mpi_calc` run, there are some processes acting as masters to distribute calculation tasks, while other processes are workers. 
One worker (called leader) plays a role of an agent between masters and other workers, while other workers launch Fortran programs to perform computationally expensive tasks. 
The number of master processes can be limited to avoid waste of resources, and one can also consider [oversubscribing](https://www.open-mpi.org/faq/?category=running#oversubscribing) for better efficiency since the masters and the leader don't perform long calculations. 

`pynfam_fit_wrapper` can restart from existing solution files, so running out of walltime limit is not a big issue. 
`pynfam_fit_warpper` records its input parameters in files (pynfam_inputs_*category*.pkl/.json), and when it's restarted it will check what parameters are different to determine whether existing hfb / fam solutions should be recalculated. 
For example, when a time-odd Skyrme parameter or isoscalar pairing strength is changed, hfb solutions will not be recalculated but fam solutions will be, because the hfb part does not depend on these parameters but the fam part does. 
Another example is that all the existing solutions will be reused if we restart with exactly the same set of input parameters. 
Therefore, it is not recommeded to delete existing solution files before we start a new function call of `pynfam_fit_warpper`, unless we encounter a bug that may destroy some old files during the prefious run. 

## Summary of PNFAM calculation

The mathematical essense of PNFAM is described below in brief, which may provide some insights for building an emulator. 
PNFAM solves a linear equation $$(A+\omega B)\vec{x} = \vec{b}(F)$$ in an iterative way for a number of $\omega$ values.  
Matrix $A$ depends linearly on parameters we want to fit, namely $$A=A_0+\sum_i c_i P_i$$ where $c_i$ is one of the parameters and matrix $P_i$ describes how $A$ depends on $c_i$. 
Matrix $B$ does not depend on any parameters. 
Vector $\vec{b}$ depends on the transition operator $F$ while the matrices on the left-hand side do not depend on it. 
The strength function $S(\omega;F)$ is defined as the inner product of $\vec{b}$ and $\vec{x}$, i.e. $b^\dagger x$. 

For GTR, $\omega$ values are written as $\omega = \Omega + i\Gamma$ where the energy $\Omega$ and half width $\Gamma$ are real-valued. $\Gamma$ is fixed while $\Omega$ = `np.arange(energy_min, eenergy_max+de, de)` with `energy_min`, `energy_max` and `de` specified by the user. 
The GT strength distribution is the imaginary part of $$-\frac{1}{\pi}\left[S(\omega;\mathrm{GT\_K0}) + 2S(\omega;\mathrm{GT\_K1})\right].$$ The maximum of the strength distribution is treated as the GTR, and the corresponding energy $\Omega$ is the model prediction of the GTR energy. 
In PNFAM the factor of $-1/\pi$ is included in the strength function, just for simplicity.

Theoretically speaking, $S(\omega;F)$ has the form of
$$S(\omega;F)=\sum_{n}\frac{S_n(F)}{\omega-E_{n}}$$
where $S_n(F)$ and $E_n$ are real-valued. 
It can be easily verified that the imaginary part of $S(\omega;F)$ has the form of
$$-\Gamma \sum_n \frac{S_n(F)}{(\Omega - E_n)^2 + \Gamma^2}$$
for $\omega = \Omega + i\Gamma$. 

For the decay half life, $S(\omega;F)$ is calculated for a number of operators $F$ and for $\omega$ values on a circular contour in the complex plane. 
The sum of contour integrations weighted by phase factors that depend on $F$ and $\omega$ gives the decay rate $\lambda$. 
Then the half life is calculated by $$T_{1/2} = \frac{\ln(2)}{\lambda}.$$

## Installation
### Fortran
* BLAS and LAPACK libraries are neede to compile Fortran programs. It's recommended to use Intel compilers and Intel MKL for best performance. 
* Make sure that all the Makefiles have the correct flags for BLAS and LAPACK before compiling. For HFBTHO it can be done by specifying the correct value for variable `COMPUTER`. 
* Run `make` in the directory "exes/hfbtho_blocking_ba9f1a8/src" to compile HFBTHO. Copy the executable "exes/hfbtho_blocking_ba9f1a8/src/hfbtho/hfbtho_main" to the directory given by the parameter "exes_dir" of `pynfam_fit_wrapper`. 
* Run `make pnfam_main.x` in the directory "exes/pnfam/" to compile PNFAM. Copy the executable "exes/pnfam/pnfam_main.x" to the directory given by the parameter "exes_dir" of `pynfam_fit_wrapper`. 
* **Now `pynfam` can search for the Fortran executables recursively in the directory (and subdirectories) given by the parameter "exes_dir"**, so you can also just use the respository contaning all the source codes and executables as "exes_dir", without the necessity to copy the executables out to a separate directory. 
  
### Python
* Python version >= 3.5 is required for `pynfam_fit_wrapper`.
* PyNFAM depends on following packages: numpy, scipy, pandas, mpi4py, f90nml. 
* One can either use the source code to run PyNFAM, or use `python setup.py install --user` (`--user` can be removed if the administration privilege is granted) to install PyNFAM. 
* When the code is changed, first uninstall PyNFAM by `python -m pip uninstall PyNFAM` before reinstall it. 

## Parameters of `pynfam_fit_warpper`

If there is no special requirement, it will be OK to keep default parameter values, except `override_setts_fit` which contains the parameters to fit. 

* **input_data**: str or pandas.DataFrame, default "pynfam_fit_input.csv".
  
  A pandas DataFrame that contains the data for fitting, or the path of a csv file with the data. In the csv file values should be separated by ',',and comments should start with '#'. 
  
  The example file "pynfam_fit_input.csv" is given with explanations of each column. It can be easily viewed in Microsoft Excel or similar softwares. 

* **categories**: str, tuple / list of str or None, default None. 

  The category (str) or categories (list / tuple) of observables to calculate. Supported categories are "GT" (Gamow-Teller resonance), "HL" ($\beta^-$ decay half life), "SD" (spin-dipole resonance) and "GAP" (pairing gap). 
  
  If None, all the categories existing in `input_data` will be calculated in the order of their first appearances in `input_data`. 
  
  Two function calls of `pynfam_mpi_calc` (one for HFB and one for FAM + GT/SD/Beta) will be launched for each category. 

* **override_setts_fit**: dict, default {}.

  A dict to generate `override_settings` for `pynfam_mpi_calc`, which provides input parameters for all the Fortran programs. 
  Some parameters in `pynfam_inputs` (`gs_def_scan`, `beta_type` and `fam_ops`) can also be specified here. 
  For each item in this parameter, its key is a str of a variable name and its value is the variable value. 

  For the prioritization of different parameter sources for `pynfam_inputs` and `override_settings`, see [Attention](##Attention). 
  <!-- **Unspecified variables will use values given in dict `default_setts` in "pynfam_default.py". **
  Some variables are not specified in "pynfam_default.py", either; they will be determined by "pynfam/config.py" if they are not given in `override_setts_fit`.  -->

  Fitting parameters are passed via `override_setts_fit`. 
  Parameters that will not change throughout the fitting are recommended to be set via `default_setts` in "pynfam_default.py", which should be put in the same directory as the script that calls the function `pynfam_fit_wrapper` so that `default_setts` can be imported. 

  A two-level dict like `default_setts` can be passed here. A one-level dict with only innermost keys is also accepted as the innermost keys are unique. 

  Recommended fitting parameters: 
  * **g0p**: $g_0^\prime$, one of Landau parameters, default 1.6 (experimental result of $^{132}$Sn, with the effective mass of UNEDF1-HFB). 
  * **g1p**: $g_1^\prime$, one of Landau parameters, default 0. 
  * h0p: $h_0^\prime$, one of Landau parameters, default 0. 
  * **vpair_t0_scaled**: isoscalar pairing strength scaled w.r.t. (divided by) the aboslute value of the isovector pairing strength, $\frac{V_{t=0}}{V_{t=1}}$, default 0. Scaling is for a better performance of $\chi^2$ minimization. 
  * **GA**: $g_A$, weak axial-vector coupling strength, for HL only, default -1. 

* **override_setts_cat**: dict, default {}.

  A dict works like `override_setts_fit`, but this one can provide different settings for different categories. 
  It can be a two-level or three-level dict. The keys in the top level are category strings. The bottom one or two levels are the same as `override_setts_fit`. 

* **return_df**: bool, default True.
  
  Whether a pandas DataFrame will be returned or not. 
  
  If True, the DataFrame of `input_data` plus columns of results will be returned. 
  See "Explanations on "input_data" and DataFrame return" below for details. 
  <!-- The column "model" contains model predictions. The column "detail" contains a series of DataFramesthe column "detail" contains a series of DataFrames: For "GT", these DataFrames give strength distributions; for "HL", these DataFrames give contributions from different orders; for "GAP", no DataFrame is given.  -->
  
  If False, the DataFrame will be converted to a dict and then returned.  
  
* **clean_files**: bool, default False.

  Whether unnecessary output files from the previous calculation will be cleaned at the end.

* **exes_dir**: str, default "./exes". 

  Path of the directory containing HFBTHO and PNFAM executables. 

* **scratch_dir**: str, default "./scratch". 

  Path of the directory to read / write solutions. Inside it, directories named after observable categories will be created (if not exist) to store solutions of different categories. 

* **backup_option**: int, default -1. 

  If 0, no backup of solutions will be done after the calculation of one category finished. 

  If +/- 1, backup will be done after the calculation of one category finished, but all the tar files that contain binary / text files directly generated by PNFAM will not be copied to save storage. 

  If +/- 2, backup will be done after the calculation of one category finished, with all the output files copied. 

  Otherwise, an error will be raised. 

  If negative, the calculation is started from an existing backup (specified by **start_from_backup**), and there is no change in input parameters, the backup will be directed to the directory of **start_from_backup** and overwrite the old backup. 

  To reduce the number of files, directories like "000002" and "meta" inside each category folder will be tarred and named after the category. 

* **backup_max**: int, default -1. 
  
  If positive (>0), old backups will be deleted if the number of backups exceed `backup_max`. 

* **backup_dir**: str, default "./backup". 

  Path of the directory containing backups. All the subdirectories will be named by integers starting from 0. Newer backups have larger subdirectory names. When old backups are deleted, their names will not be reused for new backups. 

  <!-- If an empty string is passed, `backup_option` and `start_from_backup` will be set to 0 and -1, respectively.  -->

* **start_from_backup**: int, default -1. 
 
  If not negative (>=0), the backup subdirectory named after `start_from_backup` will be copied to scratch directory and the calculation will restart from this backup. 
  <!-- It does not make sense to use this feature if `use_fam_storage` is not -1 (PNFAM is not restarting from a binary file) and `backup_option` is not 2 (binary files are not stored in backup).  -->

* **nr_parallel_calcs**: int, None, or dict with keys of category strings and values of int or None, default None.
  
  Specify the number of master processes. It must be a positive integer.
  
  <!-- If a list / tuple is given, it must have the same length as `categories` or the category list obtained from `input_data`, in order to specify the number of master processes for each category. For one category, if `nr_parallel_calc` is None, the number of master processes will be the same as the number of data points.  -->

  <!-- If a dict is given, it will be transformed to a tuple based on `categories` or the category list obtained from `input_data`.  -->

  If a dict is given, it should have keys of category strings and values of positive int or None. 
  The nr_parallel_calcs of one category is determined by the corresponding value in the dict. 
  The category not given as a key in the dict will use None as its nr_parallel_calcs.
  If duplicate category is found in the dict (possible because of the use lower and upper cases in the category string), the last one will be used. 

* **rerun_mode**: int, str, None, or dict with keys of category strings and values of int, str or None, default None.

  If None, how the program should be restarted will be determined by checking whether current parameters match parameters of the existing solution (if exists). 
  <!-- Changing parameters that will not impact the results will not make the calculation restart.  -->

  If 0, restart from the existing solution without parameter checking. 

  If 1, the same as 0 but no HFBTHO calculation will be done, without parameter checking. Not recommended for fitting. 

  If 'FAM', the same as 0 but all the PNFAM calculations will be redone, without parameter checking. 

  If 'HFB', all the HFBTHO and PNFAM calculations will be redone, without parameter checking. 

  If 'HFB_NORESTART', all the HFBTHO and PNFAM calculations will be redone, without parameter checking, but forcing both HFB and FAM programs to restart from scratch no matter what is specified in "pynfam_default.py" and "pnfam/config.py". 

  <!-- If a list / tuple containing above possibile values is given, it must have the same length as `categories` or the category list obtained from `input_data`, in order to specify the parameter for each category.  -->

  <!-- If a dict is given, it will be transformed to a tuple based on `categories` or the category list obtained from `input_data`.  -->

  If a dict is given, it should have keys of category strings and values of all the cases give above. 
  The rerun_mode of one category is determined by the corresponding value in the dict as described above. 
  The category not given as a key in the dict will use None as its rerun_mode. 
  If duplicate category is found in the dict (possible because of the use lower and upper cases in the category string), the last one will be used. 

  Otherwise, an error will be raised. 

* **ignore_hfb_nonconv**: int, None, or dict with keys of category strings and values of int or None, default 2. 

  Parameter in `pynfam_inputs['hfb_mode']['ignore_nonconv']`. Quote "README.md" below (with the case of None added). 

  > Define how to handle non-converged (nc) hfb solutions.
  >
  >  * 0 - stop for any nc soln
  >  * 1 - discard odd nc solns (stops if ANY even or ALL odd are nc)
  >  * 2 or None - discard all nc solns (stops if ALL even or ALL odd are nc)
  >  * 3 - use nc solns

  If a dict is given, it should have keys of category strings and values of all the cases give above. 
  The ignore_hfb_nonconv of one category is determined by the corresponding value in the dict as described above. 
  The category not given as a key in the dict will use None as its ignore_hfb_nonconv. 
  If duplicate category is found in the dict (possible because of the use lower and upper cases in the category string), the last one will be used. 
  
* **hfb_only**: bool, or dict with keys of category strings and values of bool, default False.

  Whether only HFB calculations will be performed. 

  If a dict is given, it should have keys of category strings and values of bool. 
  The hfb_only of one category is determined by the corresponding value in the dict as described above. 
  The category not given as a key in the dict will use None as its hfb_only. 
  If duplicate category is found in the dict (possible because of the use lower and upper cases in the category string), the last one will be used. 

<!-- * **use_fam_storage**: int, None or list/tuple of int and None, default None. 

  **To avoid confusion, this option is hidden and hard coded as None.**

  Specify how PNFAM handles binary files that containing solutions. 

  If None, PNFAM will use the value specified in "pynfam_default.py" (or "pnfam/config.py" if not specified in "pynfam_default.py"). 
  
  If 0, PNFAM will not read or write any binary file. 
  
  If 1, PNFAM will read the binary file and restart from it (if exists), but will not write to the binary file. 
  
  If -1, PNFAM will read the binary file and restart from it (if exists) and write the new solution to it; otherwise, an error will be raised.  

  If a list / tuple containing above possibile values is given, it must have the same length as `categories` or the category list obtained from `input_data`, in order to specify the parameter for each category. 

  If the binary file containing a solution that was obtained from a set of parameters close to the current parameters, setting `use_fam_storage` as -1 will make PNFAM converge faster. Otherwise, it's better to start PNFAM from scratch. 
  **Tests show that the convergence acceleration obtained by restarting from an existing solution is not significant, so it's recommended to start from scratch for all the cases.** -->

* **use_ratinterp**: int, default 1. 

  Whether a rational interpolation is used to find the GTR or SDR. This feature allows an accurate extraction of peak position with a small number of points calculated. 

  If 0, the peak position is determined by finding the maximum strength in all the computed points. 

  If 1 or 2, the inverse of the complex strength function $S(\omega)$ will be interpolated by thiele's interpolator adapted from class `phaseSpace`, with the assumption that $S(\omega^*)=S^*(\omega)$ (all the poles are on the real axis with real-valued residuals); 
  then a SciPy mimization function is called to find the peak corresponding to the GTR or SDR. 
  The Nelder-Mead method is utilized in the minization function because of its robustness; other methods may fail and give an unstable result. 

  If 2, the numerical error is also calculated by error propagation. The source of the numerical error is the FAM convergence accuracy (si) at each point, and derivatives are computed by the finite difference method. 

  Otherwise, an error will be raised. 

* **assert_consistency**: bool, default False. 
  
  Whether parameter consistency between different processes will be checked.

* **comm**: MPI communicator or None, default None. 

  Pass a MPI communicator to use. If None, MPI_COMM_WORLD (all the processes available to use) will be used. 

* **check**: bool, default False. 

  If True, no actual calculation will be done. The program will only check whether input parameters are valid, and return a boolean indicating whether all the rerun modes are 0. 
  This return can be used to quickly check whether the backup you start from has the same set of parameters as what you want. 

* ****kwargs**:

  Parameters provided here will play the same role as those given in `override_setts_fit`. 
  
  For the use of "kwargs" in Python please read https://book.pythontips.com/en/latest/args_and_kwargs.html

## Return of `pynfam_fit_wrapper`

A pandas DataFrame or a dict. The dict return is directly transformed from the DataFrame. 
<!-- The residual is the difference between the model prediction and the experimental value (if provided); the normalized residual is obtained by dividing the residual by the normalization factor "sigma" (if provided).  -->

The DataFrame containing all the results will also be stored on dist at "*scratch_dir*/result.pkl" (serialized by pickle) and "*scratch_dir*/result.csv" (without column "detail"). 
The DataFrame for each category can also be found at their own folders, with file names "result_*category*.pkl" and "result_*category*.csv". 

## Class `pynfam_residual`

This is a wrapper class to provide residuals for fitting / least_squares routines like POUNDerS. 
### Attributes:
* **args**: dict

  Parameters to be passed to `pynfam_fit_wrapper`.

* **nonConv_penalty**: float or None, default None. 

  The penalty $\lambda$ for the non-convergent result. The residual of the non-convergent case will be multiplied by a factor of $(1+\left|\lambda\right|)$. 
  If None, no penalty will be applied. 

* **comm**: MPI communicator. 
* **rank**: MPI process rank. 
* **count**: Number of method `fun` evaluations. It will start from `start_from_backup`.
* **result**: Result returned by `pynfam_fit_wrapper`, default None. 
  
* **restart_in_seq**: bool, default False. 
  Whether enable the restart of calculations from existing backups. 
  If True, the method `fun` will try to recover results from backup to avoid expensive recalculations, based on `count`. 

* **special_params_names**: 1D array of str, default (). 
  Names of special parameters (see method `set_special_params` and `fun_special`). 
* **special_params_vals**: 1D array of float, default (). 
  Current values of special parameters (see method `set_special_params` and `fun_special`). 
* **special_params_bounds**: 2-tuple of 1D arrays, default (-np.inf, np.inf). 
  Bounds of special parameters (see method `set_special_params` and `fun_special`). 
* **special_params_reeval**: function, default `lambda *args: None`. 
  Function to reevaluate the result with new values of special parameters (see method `set_special_params`, `fun_special` and `vary_GA_only`). 
* **special_params_xtol**: float, default 1e-12. 
  Tolerance for the optimization over speical parameters (see method `set_special_params` and `fun_special`). 
* **special_params_outfile**: str or None, default None. 
  File to log the optimization over special parameters (see method `set_special_params` and `fun_special`). 
  No file output if None. 

### Methods:
* **\_\_init\_\_**(self, input_data, categories, override_setts_fit, return_df, clean_files,
          exes_dir, scratch_dir, backup_dir, backup_option, backup_max, start_from_backup,
          nr_parallel_calcs, rerun_mode, ignore_hfb_nonconv, use_ratinterp, assert_consistency,
          comm, check, 
          special_params_names, special_params_x0, special_params_reeval, special_params_bounds, 
          special_params_xtol, special_params_outfile, nonConv_penalty, restart_fit_routine,
          **kwargs):

  Initilize attributes. 
  
  `self.args` takes all the input parameters except `nonConv_penalty`, `special_params_*` or `restart_fit_routine`, and store them in a dict. 
  `**self.args` will be passed as input parameters to function `pynram_fit_wrapper`. 
  The default values of parameters in `self.args` are the same as those in `pynfam_fit_wrapper`. 
  Input parameter `hfb_only` is not allowed here; its default value (False) is passed to `pynfam_fit_wrapper`. 

  The default value of `nonConv_penalty` is None; see method **extract_residual** for details.

  The default values of `special_params_*` are
  `
  special_params_names=(), special_params_x0=(), special_params_reeval=lambda *args: None, special_params_bounds=(-np.inf, np.inf), 
  special_params_xtol=1e-12, special_params_outfile=None
  `
  ; see method **set_special_params** for details.

  The default value of `restart_fit_routine` is False; see method **fun** for details.

* **set_special_params**(self, names, x0, routine, bounds=(-np.inf, np.inf), xtol=1e-12, outfile=None):

  Set names, starting point, bounds, evaluations function, tolerance and output file for special parameters; 
  names and x0 must be an 1D arrays with the same length, and bounds must be 2-tuple of 1D arrays. 
  These input arguments are stored into attributes **special_params_names**, **special_params_vals**, **special_params_bounds**, **special_params_reeval**, **special_params_xtol** and **special_params_outfile**. 

  Special parameters can be optimized over separately; see method **fun_special** for details.
  Usually parameter **GA** is set as an special parameter. 

  The argument `routine` is a function called as `routine(self.result, **dict(zip(self.special_params_names, self.special_params_vals)))`,  
  where `result` is what is returned by `pynfam_fit_wrapper` with old values of special parameters, and `x` gives an array of new values of special parameters. 
  This function should return the new result calculated with the new values of special parameters, 
  but still in the same format as `pynfam_fit_wrapper`. 
  For **GA**, the static method **vary_GA_only** should be the `routine` here.

  It is called by **\_\_init\_\_** to set attributes about special parameters. 

* **call_pynfam_fit_wrapper**(self, **fit_vars):
  
  Call function `pynfam_fit_wrapper` with `**fit_vars`, `**self.args` and `**self.kwargs` as input parameters. 
  If the same key appears in `fit_vars` and `self.args` / `self.kwargs`, only the one in `fit_vars` will take effect. 

  It can only be called by method **call_pynfam_extract_res** if `restart_fit_routine` is True; otherwise, an error will be raised. 

* **extract_cols**(result, *cols):
  
  Extract the ndarrays corresponding to columns / keys given by `cols` from `result`, the return of `pynfam_fit_wrapper`. 
  If only there is only one item in `cols`, the corresponding ndarray will be returned; otherwise, this function will return a list of ndarrays. 

  This is a static method that can be used without creating an object. 

* **call_pynfam_extract_cols**(self, *cols, **fit_vars):

  Call function `pynfam_fit_wrapper` with `**fit_vars` and class attributes as input parameters, and then extract ndarrays corresponding to columns / keys given by `cols`. 

  It cannot be called if `restart_fit_routine` is True; otherwise, an error will be raised. 

* **extract_residual**(result, nonConv_penalty):

  Extract residuals from `result`, the return of `pynram_fit_wrapper`, with (1+`nonConv_penalty`) multiplied onto non-convergent points if `nonConv_penalty` is not None. 
  It returns "norm_res" when it exists; otherwise, "residual" is returned if it exists. It returns "model" when no residual is provided. 

  This is a static method that can be used without creating an object. 

* **call_pynfam_extract_res**(self, **fit_vars):

  Call function `pynram_fit_wrapper` with `**fit_vars` and class attributes as input parameters, and then extract residuals with `self.nonConv_penalty` applied. 

  It can only be called by method **fun** if `restart_fit_routine` is True; otherwise, an error will be raised. 

* **fun**(self, xin, var_names):

  Call method `call_pynfam_extract_res` with the values of fitting variables given in array-like `xin` with shape (n,); 
  the names of these fitting variables are passed by array-like `var_names` with shape (n,). 
  If there is only one fitting variable, `xin` and `var_names` must still be array-like with shape (1,). 

  A dict with the elements in `var_names` as keys and those in `xin` are constructed and passed to `call_pynfam_extract_res` as `fit_vars`. 

  This method is suitable for fitting / least_squares routines like POUNDerS. Several methods based on `fun` are also provided below for convenience. 

* **fun_g0p**(self, xin):

  Call method `fun` with "g0p" as the fitting variable. 
  "g0p" is passed via array-like `xin` with shape (1,).

* **fun_g0p_vpair0**(self, xin):

  Call method `fun` with "g0p" and "vpair_t0_scaled" as fitting variables. 
  "g0p" and "vpair_t0_scaled" are passed via array-like `xin` with shape (2,); the first element in `xin` is "g0p" and the second is "vpair_t0_scaled". 

* **fun_...\_...\_......**(self, xin):
  
  Other methods with names in the above form perform in the same way as method `fun_g0p_vpair0`, except that the names of fitting variables passed via array-like `xin` are given in the function name after the first underscore. 
  The sequence of these variables should be observed in `xin`. 

  It should be easy to extend this class by adding more similar methods to accomodate more fitting variables (read https://www.geeksforgeeks.org/extend-class-method-in-python/ for details). 

* **fun_special**(self, xin, var_names):

  Call method `call_pynfam_extract_res` with the values of fitting variables given in array-like `xin` with shape (n,); 
  the names of these fitting variables are passed by array-like `var_names` with shape (n,). 
  In addition, least-square optimization over special parameters are done separately by calling scipy.least_squares and **special_params_reeval**; 
  special parameters should not be given via `xin` and `var_names`. 
  In other words, the residuals returned by this method corresponds to 
  $$\min_{\mathrm{special\ params}} \chi^2(\mathrm{all\ params})$$. 

  If there is only one fitting variable, `xin` and `var_names` must still be array-like with shape (1,). 

  This method is suitable for fitting / least_squares routines like POUNDerS. 

* **static_find_minChiSqr_in_backups**(var_names, backup_path='./backup', result_filename='result'):
  
  Find the solution with minimum $\chi^2$ from backups in `backup_path`, with the names of fitting variables given by an array-like `var_names`. 
  The filename (without extension) of the result pkl file to read is passed by `result_filename`. 
  If there is only one fitting variable, `var_names` should be of shape (1,). 

  Its return is a 4-tuple consisting of the parameter vector, the residual vector, the value of minimum $\chi^2$ (2-norm of the residual vector) and the corresponding backup directory number. 

  This is a static method that can be called without creating an instance. 

* **find_minChiSqr_in_backups**(var_names, result_filename='result'):

  Find the solution with minimum $\chi^2$ from backups in `self.backup_dir`, with the names of fitting variables are given by array-like `var_names`. 
  The filename (without extension) of the result pkl file to read is passed by `result_filename`. 
  If there is only one fitting variable, `var_names` should be of shape (1,). 

  Its return is a 4-tuple consisting of the parameter vector, the residual vector, the value of minimum $\chi^2$ (2-norm of the residual vector) and the corresponding backup directory number. 

  This is the non-static version of `static_find_minChiSqr_in_backups`. 

* **vary_GA_only**(df_in, GA, inplace=False): 
   
   Vary the weak axial-vector coupling **GA** and obtain new beta-decay half lives. 

   Old result (full DataFrame or dict, with "detail" column included) returned by `pynfam_fit_wrapper` is passed via `df_in`, 
   and the new value of GA is passed via `GA`. 
   New result will be returned in the same format as `df_in`. 

   If `inplace` is True, `df_in` will be modified in place, and the method will return `None`. 

   This is a static method that can be called without creating an instance. 

## Explanations on `input_data` and DataFrame return

All the columns needed in the parameter `input_data` / DataFrame return are explained below. 
Columns with all NaN will be dropped in the return. 
An example input csv file called "pynfam_fit_input.csv" is also given. 

For simplicity, columns with all NaN will be removed when the DataFrame is returned or written to disk. 

* **label**: Label of the nucleus, for human reading only. 
* **category**: Type of observable, can be "GT", "HL", "SD" or "GAP". 
* **Z**: Proton number.
* **N**: Neutron number.
* **A**: Mass number Z+N, for human reading only. 
  
* **beta_type**: 

  Beta decay type. It can be '-', '+' or 'c', corresponding to beta minus, beta plus and electron capture. 
  
  If this column is missing or a cell is left blank, default value '-' will be used. 

* **func**: 
  
  Function to be applied on raw model predictions and experimental values. 

  * For string 'LOG' / 'LOG10' / 'LN', `np.log10` / `np.log10` / `np.log` will be used. Default value for 'HL' should be 'LOG'. 
  * For string 'RATE', function $f(x)=\ln 2 /x$ will be used. 

  If it's missing or a cell is left blank, no function will be applied. 

* **exp** and **exp_raw**: 
  
  Experimental value (true value). If column "func" exists, column "exp" will be the values after the application of "func", i.e. "func(exp_raw)". 

  If only column "exp" exists in the input, it will be moved to column "exp_raw" in the output, which causes inconsistency between the input and output. 

* **energy_min** and **energy_max**: 
  
  Search region for the GTR or SDR. Ignored for other categories. 

* **half_width**: 
  
  Imaginary part of $\omega$, for GT or SD only. Ignored for other categories. 

  For GT or SD, if it's missing or a cell is left blank, the default value in "pynfam_default.py" (or "pnfam/config.py" if not specified in "pynfam_default.py") will be used. 

* **de_hw_ratio**: 
  
  de/half_width, ratio between de and half_width. 
  "de" is the grid spacing inside the search region. 
  Recommended value is <~ 1.0. 

  For GT or SD only. Ignored for other categories. 
  For GT or SD, if it's missing or a cell is left blank, the default value in "pynfam_default.py" (or "pnfam/config.py" if not specified in "pynfam_default.py") will be used. 

* **deformation**: 
  
  Deformation scan region. Default value (-2,(-0.2,0.0,0.2)) will be used if it's missing or a cell is left blank. 

  It can be a string containing letters 'S' (spherical, $\beta_2=0.0$), 'P' (prolate, $\beta_2=0.2$) and 'O' (oblate, $\beta_2=-0.2$). 

  It can also be a string in the form of a tuple or list that specifies an array of $\beta_2$, like '(0.0,0.2)', '\[0.0,0.2\]' or '0.0,0.2'. If there is only one element, strings like '0.0', '(0.0)' or '\[0.0\]' are also acceptable. 

  It can also be a tuple in the form of (-2, (-0.2,0.0,0.2)). 

  For the spherical case ($\beta_2=0$), only operator GT_K0 (SD_K0) will be calculated for category GT (SD). 

* **op**: 
  
  Operator to use, for HL only. Ignored for other categories. 

  It can be 'ALL', 'ALLOWED', 'FORBIDDEN', '0+', 'GAMOWTELLER', 'SPINDIPOLE', '1+', '0-', '1-', '2-'. 
  For HL if it's missing or a cell is left blank, the default value is "ALL". 

  If only "ALLOWED" transitions contribute to the beta decay rate, we can set op as "ALLOWED" to reduce computations.

* **sigma**:

  Normalization factor for the residual. A warning will be printed if the sigma values of the same category are different. 

  If it's missing or a cell is left blank, the corresponding "norm_res" (see below) will be left blank. 

* **model** and **model_raw**:
  
  Model prediction. If column "func" exists, column "model" will be the values after the application of "func", i.e. "func(model_raw)". 

* **num_error** and **num_error_raw**:

  Numerical error estimation of the model prediction. The source of the error is FAM convergence epsilon at each point. 

  If column "func" exists, column "model" will be the errors after the application of "func", i.e. the absolute value of func's derivative multiplied by "num_error_raw". 
  
  For GTR/SDR see the explanation at the input parameter `use_ratinterp`. 
  For HL, the error propagation is done in class `shapeFactor`. 

<!-- * **model_fun**:

  Model prediction with specified functions applied. See the input parameter `apply_fun` of function `pynfam_fit_wrapper`.

* **num_error_fun**:

  Numerical error estimation of column "model_fun", based on column "num_error" and the derivative of the applied function. 

* **exp_fun**:

  Experimental value with specified functions applied. See the input parameter `apply_fun` of function `pynfam_fit_wrapper`. -->

* **residual**:

  Difference between "model" and "exp".

* **norm_res**:

  Normalized residual, i.e. the result obtained by dividing "residual" by "sigma". 

  If column "sigma" is missing or a cell in column "sigma" is left blank, the corresponding norm_res will be left blank. 

* **convergence**:

  Whether the calculation is converged or not. 

* **detail**:

  DataFrame that contains some details about the calculation. 

## Attention

1. Prioritization of differnet sources for parameters in `pynfam_inputs` and `override_setting`:
   
   `DEFAULTS` in "pnfam/config.py" < hard-coded default values in `pynfam_fit_wrapper_root` < `default_setts` in "pynfam_default.py" < `override_setts_cat` < `override_setts_fit` < `kwargs` < `input_data`

   Please note that 
   1. `input_data` can only affect part of parameters in `pynfam_inputs` and `override_settings`. 
   2. only part of parameters have hard-coded default values.

2. Make sure that you have "pynfam_default.py" in the same directory with your script that calls `pynfam_fit_wrapper`. Otherwise, an warning will be given. 
3. One FAM binary file can take 6~8 MB. If you choose `use_fam_storage`=-1 in "pynfam_default.py" and `backup_option`=2, you may run out of storage soon. 
4. For category "GT", "SD" or "GAP", `fam_ops` specified at any place (`default_setts` in "pynfam_default.py", `override_setts_cat` or `override_setts_fit`) will have no effect. 
5. For category "GT" or "SD": 
   1. If the search region specified by "energy_min" and "energy_max" is not large enough, the resonance peak may fall out of it. 
   However, it's not recommended to use a very large region so that the computational cost is under control. 
   1. Both "energy_min" and "energy_max" are given as the excitation energy relative to the ground state of the daughter nucleus. 
   It is converted to the QRPA energy using the mass excess given in the atomic mass evaluation 2020. 
   See the appendix of Phys. Rev. C 101, 044305 (2020) for equations. 
   1. If the grid spacing "de" is not small enough, a small variation of fitting parameters may not give a different result. 
   I'm not sure whether this will bring any difficulty to the statistics. It should be noted that the error of experimental values is >~ 0.1 MeV, so from the physics perspective it is not necessary to use a very small spacing. 

## Example

"pynfam_fit_test.py" provides a script for a quick test. In this example "nr_points" is changed to 6 to reduce the computational cost for a test.
More examples are provided in the directory "examples". 

A wrapper for `pynfam_fit_wrapper` is needed if a statistics module is going to be used. For example, 
```
    def pynfam_fit_stat_wrapper(g0p):
        get MPI rank
        if rank == 0:
            prepare input parameters for pynfam_fit_wrapper
            call comm.bcast to braodcast input parameters
        else:
            call comm.bcast to receive input parameters
            if receive a termination flag: return None
        result = pynfam_fit_wrapper(...)
        post processing of result and return
    
    if __name__ == 'main':
        ...
        get MPI rank
        if rank == 0:
            solve(pynfam_fit_stat_wrapper) # call function provided by statistics package
            call comm.bast to braodcast the termination flag
        else:
            while True:
                result = pynfam_fit_stat_wrapper(None)
                if result is None: break
        ...
```
See "fit_examples/pounders" for the use of `pynfam_residual` with POUNDerS. 
For the use of POUNDerS, read comments in "fit_examples/pouders.py". The Python wrapper of POUNDerS makes the use of POUNDerS similar to [`scipy.optimize.least_squares`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.least_squares.html). 

<!-- It's also good to use MPI Spawn function in `pynfam_fit_stat_wrapper` for dynamic process management. Pseudocode remains to be provided.  -->
