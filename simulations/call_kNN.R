# Call kNN
library(reticulate)

# Use appropriate Python environment
use_virtualenv("r-reticulate", required = FALSE)  # or use_python("/usr/bin/python3")

# Ensure required Python packages are available
py_install(c("scikit-learn"), pip = TRUE)

call_kNN_dist <- function(D, y, k, precomputed=1){

np <- import("numpy")

#D <- r_to_py(D, convert = FALSE)
#py$D <- D
#py$y <- r_to_py(as.integer(y))
#py$k <- r_to_py(as.integer(k))

py$D <- np$array(D)
py$y <- np$array((as.integer(y)))
py$k <- as.integer(k)
py$precomputed <- as.integer(precomputed)

#py$k_max <- r_to_py(15)

# Call sklearn
py_run_string("

import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.neighbors import KNeighborsClassifier


# indices only
n = D.shape[0]
idx = np.arange(n)

train_idx, test_idx = train_test_split(
    idx,
    test_size=0.2,
    #stratify=y
)

if precomputed == 1:
    # slice distance matrix correctly
    D_train = D[np.ix_(train_idx, train_idx)]
    D_test  = D[np.ix_(test_idx, train_idx)]

if precomputed == 0:
    D_train = D[train_idx, ]
    D_test  = D[test_idx, ]


y_train = y[train_idx]
y_test  = y[test_idx]

# KNN with precomputed distances
if precomputed == 1:
    knn = KNeighborsClassifier(
        n_neighbors=k,
        metric='precomputed'#,
        #weights='distance'
    )
# KNN with precomputed distances
if precomputed == 0:
    knn = KNeighborsClassifier(
        n_neighbors=k#,
        #weights='distance'   
    )

knn.fit(D_train, y_train)
y_pred = knn.predict(D_test)
")

 # bring results back to R
res = list(
    y_test = py$y_test,
    y_pred = py$y_pred
  )

return(res)

}

# The call would be:
#res <- call_kNN_dist(D, y, k = 5)

#mean(res$y_test == res$y_pred)