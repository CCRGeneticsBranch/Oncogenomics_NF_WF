#!/usr/bin/env python3

import sys

def check_strandedness(filename, tool):
    with open(filename, 'r') as f:
        header = f.readline().strip().split('\t')
        line = f.readline().strip().split('\t')
        forward_pct = float(line[header.index('ForwardPct')])

        if forward_pct > 0.75:
            if tool == 'picard':
                return 'FIRST_READ_TRANSCRIPTION_STRAND'
            else:
                return 'forward'
        elif forward_pct < 0.25:
            if tool == 'picard':
                return 'SECOND_READ_TRANSCRIPTION_STRAND'
            else:
                return 'reverse'
        else:
            if tool == 'picard':
                return 'NONE'
            else:
                return 'none'

if __name__ == '__main__':
    filename = sys.argv[1]
    tool = sys.argv[2] if len(sys.argv) > 2 else 'rsem'
    result = check_strandedness(filename, tool)
    print(result)

