BEGIN {
    first = 0
}

!/^#/ && NF > 0 {
    if (first == 0) {
        max = $NF
        min = $NF
        first = 1
        argmax = $1
        argmin = $1
    }
    if (first > 0) {
        if (max < $NF) {
            max = $NF
            argmax = $1
        }
        if (min > $NF) {
            min = $NF
            argmin = $1
        }
    }
}

END {
    printf("%.7f\n", argmin)
}
