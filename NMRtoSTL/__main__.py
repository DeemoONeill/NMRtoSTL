try:
    from NMRtoSTL.NMRtoSTL import main
except ModuleNotFoundError:
    from NMRtoSTL import main
import click


@click.command()
@click.argument("filename", nargs=1)
@click.option(
    "-o",
    "--outfile",
    default=None,
    type=str,
    help="filename and location for the output file",
)
@click.option(
    "--f1",
    nargs=2,
    type=float,
    help="min and max values for f1 passed in together",
)
@click.option(
    "--f2",
    nargs=2,
    type=float,
    help="min and max values for f2 passed in together",
)
@click.option(
    "--f1min",
    default=None,
    type=float,
    show_default=True,
    help="Minimum value on the f1 axis",
)
@click.option(
    "--f1max",
    default=None,
    type=float,
    show_default=True,
    help="Maximum value on the f1 axis",
)
@click.option(
    "--f2min",
    default=None,
    type=float,
    show_default=True,
    help="Minimum value on the f2 axis",
)
@click.option(
    "--f2max",
    default=None,
    type=float,
    show_default=True,
    help="Maximum value on the f2 axis",
)
@click.option(
    "-s",
    "--stack",
    default=-1,
    type=int,
    help="Number of 1D spectra to stack together",
)
@click.option(
    "-S",
    "--sigma",
    default=(1, 10),
    nargs=2,
    type=int,
    help="""sigma : pair of ints, optional
            Sigma values for x and y axis to use for gaussian smoothing of peaks.
            The default is (0,0)""",
)
@click.option(
    "-z",
    "--size",
    default=(10, 10, 6),
    nargs=3,
    type=float,
    help="""X, Y and Z dimensions of final mesh (in mm).
            The default is (5, 5, 5).""",
)
@click.option(
    "-t",
    "--thickness",
    default=0.5,
    nargs=1,
    type=float,
    help="""Thickness of the base (in mm)""",
)
@click.option(
    "-T",
    "--threshold",
    default=0,
    nargs=1,
    type=float,
    help="""Threshold to remove noise from baseline as a percentage of the
            height of the tallest peak. The default is 0%.""",
)
@click.option("-v", "--verbose")
def collect_args(
    filename,
    f1,
    f2,
    f1min,
    f1max,
    f2min,
    f2max,
    stack,
    outfile,
    sigma,
    size,
    thickness,
    threshold,
    verbose,
):
    if f1:
        f1min, f1max = f1
    if f2:
        f2min, f2max = f2
    process_args = dict(sigma=sigma, size=size, threshold=threshold)
    main(
        filename,
        f1min,
        f1max,
        f2min,
        f2max,
        stack,
        outfile,
        thickness,
        process_args,
        verbose,
    )


collect_args()
