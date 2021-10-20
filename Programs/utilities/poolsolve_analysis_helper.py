# helper functions for analyzing the reasons why there are multiple pool entries of the same dimension for many rows
import numpy as np
import matplotlib.pyplot as plt


# analyzes if misreadings or mutations are to blame
def examine_mutations(all_combos, barcodeMap):

    # counts everything
    tots_by_type = {}  # key is tuple of pool types, value is array of number of occurrences by hamming distance
    counts_by_type = {}

    types = ['r', 'c', 'pr', 'pc']
    for i in range(len(types)):
        for j in range(i, len(types)):
            tots_by_type[(types[i], types[j])] = np.zeros(9)
            counts_by_type[(types[i], types[j])] = np.zeros(9)

    for combo in all_combos:
        type1 = get_pool_type(combo[0])
        type2 = get_pool_type(combo[1])
        ham = get_ham(barcodeMap[combo[0]], barcodeMap[combo[1]])
        tup = (type1, type2)
        if tup not in tots_by_type:
            tup = (type2, type1)
        tots_by_type[tup][ham] += all_combos[combo]
        counts_by_type[tup][ham] += 1

    means_by_type = {}
    for tup in tots_by_type:
        for i in range(len(counts_by_type[tup])):
            counts_by_type[tup][i] = max(1, counts_by_type[tup][i])
        means_by_type[tup] = tots_by_type[tup] / counts_by_type[tup]

    # plots counts as function of hamming distances for each pool type
    for type1 in types:
        plt.figure()
        leg = []
        for type2 in types:
            tup = (type1, type2)
            if tup not in means_by_type:
                tup = (type2, type1)
            plt.plot(range(9), means_by_type[tup])
            leg.append(type1 + '_' + type2)
        plt.legend(leg)
        plt.title('plots for ' + type1)
        plt.xlabel('hamming distance')
        plt.ylabel('count')
    plt.show()
    return


# analyzes if multiple pickings of the same mutation from a mating is to blame
def examine_mutliple_pickings(poolPresenceTable, sudokuGridLookupDict, prPools, pcPools, sudokuPoolColumns, thresh):
    num_same = 0
    num_different = 0

    num_mult_rows = 0
    num_mult_cols = 0
    num_mult_prows = 0
    num_mult_pcols = 0

    num_no_rows = 0
    num_no_cols = 0
    num_no_prows = 0
    num_no_pcols = 0

    for row in poolPresenceTable:
        prows = []
        pcols = []
        rows = []
        cols = []
        for i in range(1, len(row)):
            if row[i] >= thresh:
                if sudokuPoolColumns[i] in prPools:
                    prows.append(sudokuPoolColumns[i])
                elif sudokuPoolColumns[i] in pcPools:
                    pcols.append(sudokuPoolColumns[i])
                elif get_pool_type(sudokuPoolColumns[i]) == 'r':
                    rows.append(sudokuPoolColumns[i])
                else:
                    cols.append(sudokuPoolColumns[i])

        if len(rows) == 0:
            num_no_rows += 1
        elif len(rows) > 1:
            num_mult_rows += 1

        if len(cols) == 0:
            num_no_cols += 1
        elif len(cols) > 1:
            num_mult_cols += 1

        if len(prows) == 0:
            num_no_prows += 1
        elif len(prows) > 1:
            num_mult_prows += 1

        if len(pcols) == 0:
            num_no_pcols += 1
        elif len(pcols) > 1:
            num_mult_pcols += 1

        if len(prows) == 0 or len(pcols) == 0:
            continue

        if len(prows) == 1 and len(pcols) == 1:
            continue

        mating = sudokuGridLookupDict[prows[0]][pcols[0]].mating
        pools_seen = set()
        for prow in prows:
            for pcol in pcols:
                if sudokuGridLookupDict[prow][pcol].mating == mating:
                    pools_seen.add(prow)
                    pools_seen.add(pcol)
        if len(pools_seen) == len(prows) + len(pcols):
            num_same += 1
        else:
            num_different += 1

    print(('num_same', 'num_different'))
    print((num_same, num_different))
    print(('no rows', 'no cols', 'no prows', 'no pcols'))
    print((num_no_rows, num_no_cols, num_no_prows, num_no_pcols))
    print(('mult rows', 'mult cols', 'mult prows', 'mult pcols'))
    print((num_mult_rows, num_mult_cols, num_mult_prows, num_mult_pcols))
    return


# analyzes if contamination between pools is to blame
def examine_contamination(all_combos, individual_counts):

    # plots individual counts
    row_counts = []
    col_counts = []
    pr_counts = []
    pc_counts = []
    for type in individual_counts:
        if get_pool_type(type) == 'r':
            row_counts.append(individual_counts[type])
        elif get_pool_type(type) == 'c':
            col_counts.append(individual_counts[type])
        elif get_pool_type(type) == 'pr':
            pr_counts.append(individual_counts[type])
        else:
            pc_counts.append(individual_counts[type])

    plt.figure()
    plt.hist(row_counts)
    plt.title('row counts')

    plt.figure()
    plt.hist(col_counts)
    plt.title('col counts')

    plt.figure()
    plt.hist(pr_counts)
    plt.title('prow counts')

    plt.figure()
    plt.hist(pc_counts)
    plt.title('pcol counts')

    plt.show()

    listed_combos = []
    for combo in all_combos:
        listed_combos.append((combo, all_combos[combo], individual_counts[combo[0]], individual_counts[combo[1]]))

    listed_combos.sort(key=lambda comb: comb[1], reverse=True)
    print(listed_combos)

    return


# helper functions for everything else


# gets 'types' (col, row, prow, pcol) for a given string
def get_pool_type(pool):
    if pool.isdigit():
        return 'c'
    elif pool.isalpha():
        return 'r'
    else:
        return pool[:2].lower()


# gets hamming distance
def get_ham(code1, code2):
    ham = 0
    for c1, c2 in zip(code1, code2):
        if c1 != c2:
            ham += 1
    return ham
