import pysam

def append_duplex_tag_read_name(infile, bamindex, outfile, threads):
    if infile.endswith('.bam'):
        in_mode = 'rb'
    else:
        in_mode = 'r'
    if outfile.endswith('.bam'):
        out_mode = 'wb'
    else:
        out_mode = 'w'

    # use header from prev. file as template
    input = pysam.AlignmentFile(infile, in_mode,
                  check_header=False,
                  check_sq=False,
                  index_filename= bamindex)
    output = pysam.AlignmentFile(outfile, out_mode, template=input)

    for read in input.fetch():
        tag = read.get_tag('dx')
        new_qname = read.query_name+';'+str(tag)
        read.query_name = new_qname
        output.write(read)

# def get_duplex_tag_df(fname, threads):
#     """
#     Get a table that has the read ID and duplex tag
#     for each read in an alignment file
#     """
#     if fname.endswith('.bam'):
#         in_mode = 'rb'
#     else:
#         in_mode = 'r'
#     input =  pysam.AlignmentFile(fname, in_mode, threads=threads)
#     duplex_tag = []
#     read_ids = []
#     for read in input:
#         duplex_tags.append(read.get_tag('dx'))
#         read_ids.append(read.query_name)
#     df = pd.DataFrame()
#     df['read_id'] = read_ids
#     df['duplex_tag'] = duplex_tags
#     return df
#
#     # filter based on duplex status
#
#
#     return df






    # def get_duplex_tag_df(infile, outfile, threads):
    #     """
    #     Append the simplex tag to the bam file
    #     """
    #
    #     if infile.endswith('.bam'):
    #         in_mode = 'rb'
    #     else:
    #         in_mode = 'r'
    #     if outfile.endswith('.bam'):
    #         out_mode = 'wb'
    #     else:
    #         out_mode = 'w'
    #
    #     # use header from prev. file as template
    #     input = pysam.AlignmentFile(infile, in_mode)
    #     output = pysam.AlignmentFile(outfile, outmode, template=input)
    #
    #     for read in input:
    #         tag = read.get_tag('dx'))
    #         new_qname = read.query_name+';'+tag
    #         read.query_name = new_qname
    #         output.write(read)
