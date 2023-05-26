import sys
import argparse

from netZooPy.puma import Puma
import pandas as pd


def run_puma(args):
    '''
    Runs PUMA object creation and exports output to file.
    '''
    # Load the data as dataframes
    exprs_df = pd.read_csv(args.exprs, index_col=0, header=0, sep="\t")
    motif_df = pd.read_csv(args.motif, header=None, sep="\t")

    # Adding headers for the PUMA obj to read
    motif_df.columns = ['source', 'target', 'weight']

    # subset the expression dataframe to retain only the top NMAX
    # by mean expression. Otherwise, memory consumption is too much.

    # covering a very fringe case here where this column might
    # already be in the matrix. Just keep adding underscores to
    # create a unique column name for the row-mean values.
    mean_col_name = '__mean__'
    while mean_col_name in exprs_df.columns:
        mean_col_name = '_' + mean_col_name + '_'
    exprs_df[mean_col_name] = exprs_df.apply(lambda x: x.mean(), axis=1)

    # retain only the top NMAX and drop that mean value column since we're done with it.
    exprs_df = exprs_df.nlargest(args.nmax, mean_col_name)
    exprs_df.drop(mean_col_name, axis=1, inplace=True)

    # we need to provide a file or list of the miRna (for historical reasons with PANDA)
    # We simply need to grab the first column of the motif dataframe since that dataframe
    # does not contain any 'regular' TFs
    mirna_list = motif_df['source'].unique().tolist()

    # save_memory = False results in outputting the PUMA network in edge format
    # save_memory = True results in a matrix format
    # Pass the pandas dataframes directly rather than the PATHs.
    puma_obj = Puma(
        exprs_df,
        motif_df,
        # if ppi is None- this will set ppi to an identity matrix:
        None,
        mirna_list,
        keep_expression_matrix=True,
        save_tmp=True,
        save_memory=False
    )
    # Pull PUMA network out of object
    out_mtx = pd.DataFrame(puma_obj.puma_network)

    # Get motif order
    motif_names_ordered = motif_df['source'].drop_duplicates(keep="first")
    
    # Set headers and rownames to PUMA network
    out_mtx.set_index(motif_names_ordered, inplace=True)
    out_mtx.columns = exprs_df.index.values
    out_mtx.transpose().to_csv(
        args.output,
        sep="\t",
        header=True,
        index=True
    )


def main():
    parser = argparse.ArgumentParser(
        description="Runs PUMA on input count matrix."
    )

    parser.add_argument(
        "--motif",
        metavar="TSV",
        required=False,
        help="Motif data"
    )

    parser.add_argument(
        "--nmax",
        type=int,
        metavar="INT",
        required=False,
        default=25000,
        help="Max number of genes to consider. We take the top N by row-mean."
    )

    parser.add_argument(
        '--output',
        required=True,
        help='The name of the output file.'
    )

    parser.add_argument(
        "exprs",
        metavar="TSV",
        help="Expression count matrix"
    )

    args = parser.parse_args()

    # Run PUMA
    run_puma(args)


if __name__ == "__main__":
    main()
