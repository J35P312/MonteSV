# MonteSV
monte carlo simulation of SV breakpoint distance and gene content

1: clone the repository

2: cd MonteSV

3: python bootstrap_fusion.py > output.txt


the output file will contain every simulated rearangement on the following format:


A E C D F B G
+ + - + + - +
A: none_ "POU5F1" 1-31143120
E:  "SMAP1"_none 71528323-114054170
C:  "ITPR3"_none 33664160-58629165
D: none_ "SMAP1" 58629166-71528322
F: none_none 114054171-133382344
B:  "POU5F1"_ "ITPR3" 31143121-33664159
G: none_none 133382345-171115067

A-E: "POU5F1"- "SMAP1":
17039870.666666668


the first line gives the fragment order

the second line the orientation (+ normal, - inverted)

line 7 to 8 indicate the fragment position, and if the breakpoints are located within genes

Note that the script only simulate balanced rearrangements.

line 10 and onward presents the detected fusions (only one in this case).

the last line indicate the average distance between breakpoints.


The last line presents the p values, you may reach these by scrolling through the file, or by using tail:

tail -n 10 output.txt

