import sys

inf = open(sys.argv[1], 'r', encoding='big5-hkscs')
out = open(sys.argv[2], 'w', encoding='big5-hkscs')

table = dict()
for line in inf:
    s = line.split(' ')
    character = s[0]
    zhuyin = s[1].split('/')
    for x in zhuyin:
        if x[0] in table:
            if character not in table[x[0]]:
                table[x[0]].append(character)
            continue
        table[x[0]] = [character]
    table[character] = [character]

for key in sorted(table.keys()):
    out.write(key + '\t' + ' '.join(table[key]) + '\n')
inf.close()
out.close()