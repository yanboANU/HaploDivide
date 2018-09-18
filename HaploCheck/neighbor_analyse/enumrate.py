import os
import sys
import stat

a = ['A', 'T', 'C', 'G']
ans = []
for f1 in a:
    temp = ""
    temp = temp + f1
    for f2 in a:
        temp = temp + f2
        for f3 in a:
            temp = temp + f3
            for f4 in a:
                temp = temp + f4
                for f5 in a:
                    temp = temp + f5
                    ans.append(temp)
print len(ans)

stat.stat_pattern(ans)
