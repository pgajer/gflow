# X coercing 

- if (!Rf_isReal(s_X))
-     Rf_error("s_X must be numeric real.");
- PROTECT(s_X = Rf_coerceVector(s_X, REALSXP)); nprot++;
+ if (TYPEOF(s_X) != REALSXP)
+     Rf_error("s_X must be a double (REALSXP) matrix.");
+ PROTECT(s_X); nprot++;  // no copy, just protect

