import os


code = open('./keygen_original.c', 'r').read()
new_code = code

index1 = code.index('if (!compute_w(LOGN, f, g, tF, tG, tw, ')
index2 = code.index('memcpy(F, tF, DIM_N);')
part1 = code[:index1]
part2 = code[index2:]

compile_command = 'make -j16 && ./main'
# past: 41, 22, 30, 6, 10
# our:  37, 20, 25, 3, 12
limit_list = [[37, 45], [18, 26], [25, 35], [3, 10], [5, 15]]

index = 0
for i1 in range(limit_list[0][0], limit_list[0][1]+1):
    for i2 in range(limit_list[1][0], limit_list[1][1]+1):
        for i3 in range(limit_list[2][0], limit_list[2][1]+1):
            for i4 in range(limit_list[3][0], limit_list[3][1]+1):
                for i5 in range(limit_list[4][0], limit_list[4][1]+1):
                    index += 1
                    tmp_list = [i1, i2, i3, i4, i5]
                    new_code = part1 + 'if (!compute_w(LOGN, f, g, tF, tG, tw, ' + ','.join([str(vi) for vi in tmp_list]) + '))\n{continue;}\n' + part2
                    open('keygen.c', 'w').write(new_code)
                    tmp = os.popen(compile_command)
                    print(index)
                    print(tmp_list)
                    right = False
                    all_out = ''
                    for i in range(100):
                        tmp_out = tmp.readline()
                        all_out += tmp_out + '\n'
                    # if 'w = PRz([-26855,32573,-31622,24319,14895' not in all_out:
                    #     print('???')
                    #     print(all_out)
                    #     input()
                    # if '-6,60,-56,112,108,6,-55,49,-34,5,76,-104,41,-3,122,108,-11,74,-77,68,-74,-117,-117,17,-101,94,28,90,-121,49,-34,35,' not in all_out:
                    #     print('???')
                    #     print(all_out)
                    #     input()
                    first_sp_start = all_out.index('first sp:') + 10
                    first_sp_end =  all_out.index('second sp:')
                    first_sp = []
                    for vi in all_out[first_sp_start: first_sp_end].split(','):
                        try:
                            first_sp.append(int(vi))
                        except ValueError:
                            continue
                    assert len(first_sp) == 512
                    # print(first_sp)
                    if all([vi in [0, 1] for vi in first_sp]):
                        print("???")
                        print(all_out)
                        input()
# [37, 20, 25, 3, 12] wp all zero, but first sp correct