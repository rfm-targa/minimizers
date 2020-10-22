#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

Compute the set of minimizers for DNA or Protein sequences.

Code documentation
------------------
"""


import pickle
import argparse
from multiprocessing import Pool

from Bio import SeqIO


def write_to_file(text, output_file, write_mode, end_char):
    """ Writes a single text string to a file.

        Parameters
        ----------
        text : str
            Text to write to file.
        output_file : str
            Path to the output file.
        write_mode : str
            Write mode, `w` for write and `a` for append.
            `w` mode overwrites contents of files.
        end_char : str
            Char to add to end of file.
    """

    with open(output_file, write_mode) as out:
        out.write(text+end_char)


def pickle_dumper(output_file, content):
    """ Use the Pickle module to serialize an object.

        Parameters
        ----------
        content : type
            Variable that refers to the object that
            will be serialized and written to the
            output file.
        output_file : str
            Path to the output file.
    """

    with open(output_file, 'wb') as po:
        pickle.dump(content, po)


def import_sequences(fasta_path):
    """ Imports sequences from a FASTA file.

        Parameters
        ----------
        fasta_path : str
            Path to a FASTA file.

        Returns
        -------
        seqs_dict : dict
            Dictionary that has sequences ids as keys and DNA
            sequences as values.
    """

    records = SeqIO.parse(fasta_path, 'fasta')
    seqs_dict = {rec.id: str(rec.seq.upper()) for rec in records}

    return seqs_dict


def sequence_kmerizer(sequence, k_value, offset=1, position=False):
    """ Decomposes a sequence into kmers.

        Parameters
        ----------
        sequence : str
            Sequence to divide into kmers.
        k_value : int
            Value for the size k of kmers.
        offset : int
            Value to indicate offset of consecutive kmers.
        position : bool
            If the start position of the kmers in the sequence
            should be stored.

        Returns
        -------
        kmers : list
            List with the kmers determined for the input
            sequence. The list will contain strings if
            it is not specified that positions should be
            stored and tuples of kmer and start position
            if the position is stored.
    """

    if position is False:
        kmers = [sequence[i:i+k_value]
                 for i in range(0, len(sequence)-k_value+1, offset)]
    elif position is True:
        kmers = [(sequence[i:i+k_value], i)
                 for i in range(0, len(sequence)-k_value+1, offset)]

    return kmers


def generate_windows(kmers, adjacent_kmers):
    """ Create groups of adjaent kmers.

        Parameters
        ----------
        kmers : list
            List with the kmers determined for a sequence.
        adjacent_kmers : int
            Number of adjacent kmers per group.

        Returns
        -------
        windows : list
            A list with one sublist per group of adjacent
            kmers. Groups are created based on a sliding window
            whose range is incremented by 1.
    """

    window_pos = range(0, len(kmers)-adjacent_kmers+1, 1)
    windows = [kmers[i:i+adjacent_kmers] for i in window_pos]

    return windows


def determine_minimizers(seqid, sequence, adjacent_kmers, k_value,
                         position=False):
    """ Determine the minimizers for a sequence based on
        lexicographical order.

        Parameters
        ----------
        seqid : str
            Sequence identifier/header.
        sequence : str
            String representing the sequence.
        adjacent_kmers : int
            Window size value. Number of adjacent kmers per group.
        k_value : int
            Value of k for kmer size.
        position : bool
            If the start position of the kmers in the sequence
            should be stored.

        Returns
        -------
        A list with the following elements:
            seqid : str
                Sequence identifier/header.
            minimizers : list
                A list with the set of minimizers determined
                for the input sequence.
    """

    # break sequence into kmers
    kmers = sequence_kmerizer(sequence, k_value, position=position)

    # create groups of kmers for each window
    windows = generate_windows(kmers, adjacent_kmers)
    # pick smallest kmer for each window
    # min function will select minimal kmer if we pass tuples with
    # kmers and positions because it starts by comparing the kmers
    # lambda or sort are slower
    minimizers = [min(w) for w in windows]

    return [seqid, minimizers]


def default_helper(args):
    """ Helper function to pass arguments to the
        `determine_minimizers` function.

        Parameters
        ----------
        args : list
            List with list arguments to pass to function.

        Returns
        -------
        The result of the `determine_minimizers` function.
    """

    return determine_minimizers(*args)


def determine_minimizers_skipper(seqid, sequence, adjacent_kmers,
                                 k_value, position=False):
    """ Determine the minimizers for a sequence based on
        lexicographical order. Skips windows that
        cannot have a minimizer based on the minimizer
        computed in the previous iteration.

        Parameters
        ----------
        seqid : str
            Sequence identifier/header.
        sequence : str
            String representing the sequence.
        adjacent_kmers : int
            Window size value. Number of adjacent kmers per group.
        k_value : int
            Value of k for kmer size.
        position : bool
            If the start position of the kmers in the sequence
            should be stored.

        Returns
        -------
        A list with the following elements:
            seqid : str
                Sequence identifier/header.
            minimizers : list
                A list with the set of minimizers determined
                for the input sequence.
    """

    # break sequence into kmers
    kmers = sequence_kmerizer(sequence, k_value, position=position)

    i = 0
    previous = None
    sell = False
    minimizers = []
    # determine total number of windows
    last_window = (len(kmers)-adjacent_kmers)
    while i <= last_window:
        # get kmers in current window
        window = kmers[i:i+adjacent_kmers]
        # pick smallest kmer as minimizer
        minimizer = [min(window)]
        # get position in window of smallest minimizer
        minimizer_idx = window.index(minimizer[0])
        # sliding window that does not include last minimizer
        if previous is None:
            # simply store smallest minimizer
            minimizers.extend(minimizer)
        # sliding window includes last minimizer because we
        # skipped some sliding windows
        else:
            # check if minimizer is the same as the one picked
            # in the last window
            # Do not store minimizer if it is the same
            if minimizer[0] != previous:
                # get kmers smaller than last minimizer
                skipped = window[1:minimizer_idx]
                # determine if any of the smaller kmers is
                # the minimizer of a skipped window
                minimal = previous
                for m in skipped:
                    if m < minimal:
                        minimizer.append(m)
                        minimal = m
                minimizers.extend(minimizer)

        # slide by 1 if minimizer has index 0 in window
        if minimizer_idx == 0:
            i += 1
            previous = None
        # skip sliding windows based on minimizer position
        else:
            i += minimizer_idx
            # if adding minimizer index surpasses last window value we
            # might miss one last minimizer because it will fail the condition
            if i > last_window and sell is False:
                i = last_window
                sell = True
            previous = minimizer[0]

    return [seqid, minimizers]


def skipper_helper(args):
    """ Helper function to pass arguments to the
        `determine_minimizers_skipper` function.

        Parameters
        ----------
        args : list
            List with list arguments to pass to function.

        Returns
        -------
        The result of the `determine_minimizers_skipper` function.
    """

    return determine_minimizers_skipper(*args)


def map_async_parallelizer(inputs, function, threads):
    """ Create multiple processes to run a function in parallel
        for multiple inputs.

        Parameters
        ----------
        inputs : list
            List with one sublist per input.
        function
            Function that will be parallelized.
        threads : int
            Number of precesses to create.

        Returns
        -------
        results : list
            List with the returned values for each
            function call.
    """

    results = []
    pool = Pool(threads)
    process = pool.map_async(function, inputs,
                             callback=results.extend)  # chunksize=chunksize)
    process.wait()

    return results


def main(input_file, output_file, output_format, k_value,
         window_size, position, mode, threads):

    # divide big files or read per chunk to reduce memory usage
    # pass files or chunks to compute minimizers and append
    # results to output file
    # pickle.dump() can be used to save several serialized objects
    # to the same file (will have to use pickle.load() several times
    # to get all objects from file)

    # read input records
    records = import_sequences(input_file)

    # determine minimizers
    inputs = [(recid, seq, window_size, k_value, position)
              for recid, seq in records.items()]
    if mode == 'default':
        minimizers = map_async_parallelizer(inputs,
                                            default_helper,
                                            threads)
    elif mode == 'skipper':
        minimizers = map_async_parallelizer(inputs,
                                            skipper_helper,
                                            threads)

    # save results
    if output_format == 'binary':
        pickle_dumper(output_file, minimizers)
    elif output_format == 'csv':
        # create csv lines
        if position is True:
            csv_lines = ['{0},{1}'.format(m[0], ','.join(['{0}:{1}'.format(s[0], s[1]) for s in m[1]]))
                         for m in minimizers]
        elif position is False:
            csv_lines = ['{0},{1}'.format(m[0], ','.join(m[1]))
                         for m in minimizers]

        csv_text = '\n'.join(csv_lines)

        write_to_file(csv_text, output_file, 'w', '\n')


def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-i', type=str, required=True,
                        dest='input_file',
                        help='Path to a FASTA file with DNA or '
                             'protein sequences.')

    parser.add_argument('-o', type=str, required=True,
                        dest='output_file',
                        help='Path to output to save results.')

    parser.add_argument('--of', type=str, required=False,
                        choices=['binary', 'csv'], default='csv',
                        dest='output_format',
                        help='Output file format. `binary` will '
                             'use pickle to save results into a '
                             'binary file. `csv` will save results '
                             'to a CSV file with sequence identifiers '
                             'in the first field followed by a minimizer '
                             'per field.')

    parser.add_argument('--k', type=int, required=False,
                        default=4, dest='k_value',
                        help='Value for the size of the kmers.')

    parser.add_argument('--w', type=int, required=False,
                        default=4, dest='window_size',
                        help='Window size value. Minimizers will '
                             'be determined based of groups of `w` '
                             'adjacent kmers.')

    parser.add_argument('--p', action='store_false', required=False,
                        dest='position',
                        help='If the start position of the kmers '
                             'should be stored.')

    parser.add_argument('--m', type=str, required=False,
                        choices=['default', 'skipper'],
                        default='default', dest='mode',
                        help='Defines function that will be used '
                             'to determine minimizers.')

    parser.add_argument('--t', type=int, required=False,
                        default=1, dest='threads',
                        help='Number of CPU cores/threads used to '
                             'parallelize minimizer computation.')

    args = parser.parse_args()

    return [args.input_file, args.output_file, args.output_format,
            args.k_value, args.window_size, args.position, args.mode,
            args.threads]


if __name__ == '__main__':

    args = parse_arguments()

    main(args[0], args[1], args[2], args[3],
         args[4], args[5], args[6], args[7])
