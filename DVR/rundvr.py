'''script that runs a DVR'''

from .DVR import *

def rundvr(dvr_name, **pars):

    my_dvr=DVR(dvr_name, **pars)
    return my_dvr.run()


if __name__ == '__main__':
    # handle parameters passed
    import argparse

    parser = argparse.ArgumentParser(prog='rundvr')
    parser.add_argument('dvr_name',
                        type=str,
                        default='ColbertMiller1D',
                        dest='dvr_name',
                        help='the name of the DVR class to run')
    parser.add_argument('-h', dest='help_me_pls', default=False, type=bool)
    parser.add_argument('--domain', type=list, help='the range of the DVR')
    parser.add_argument('--divs', type=list, help='the divisions of the DVR')
    parser.add_argument('--pot', type=str, help='file containing the potential to use')

    args = parser.parse_args()
    if args.help_me_pls is True:
        parser.print_help()
    else:
        rundvr(args.dvr_name, **vars(args))
