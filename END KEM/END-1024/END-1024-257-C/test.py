import os
import time


code = open('keygen_original.c', 'r').read()
new_code = code

index1 = code.index('const ntru_profile SOLVE_END_257_1024 = ')
index2 = code.index('/** Compute the vector w:')
part1 = code[:index1]
part2 = code[index2:]

compile_command1 = 'make -j16 && ./main'
compile_command2 = './main'
our_data_head = "const ntru_profile SOLVE_END_257_1024 = {257, 10, 10, { 1, 1, 1, 1, 2, 4,  8, 15, 29, 57, 114}, { 1, 1, 2, 3, 6, 11, 22, 43, 84, 167  }, "
our_data_tail = ",10, { 0, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31 }, { 0, 0, 1, 2, 2, 2, 2, 2, 2, 3, 3 }};"

# { 1, 1, 2, 2, 3, 3,  4,  5,  5,  7  },
# { 1, 1, 1, 2, 3, 3,  3,  4,  5,  12 },

min_time = 100
top_words_limit = [[1, 4], [2, 8], [3, 10], [3, 10], [3, 12], [4, 14], [4, 14], [4, 16]]
# top_words_limit = [[1, 2], [2, 4], [3, 5], [3, 5], [3, 6], [4, 7], [5, 7], [12, 15]]
tmp_list = [2, 2]
for i1 in range(top_words_limit[0][0], top_words_limit[0][1]+1):
    for i2 in range(top_words_limit[1][0], top_words_limit[1][1]+1):
        for i3 in range(top_words_limit[2][0], top_words_limit[2][1]+1):
            for i4 in range(top_words_limit[3][0], top_words_limit[3][1]+1):
                for i5 in range(top_words_limit[4][0], top_words_limit[4][1]+1):
                    for i6 in range(top_words_limit[5][0], top_words_limit[5][1]+1):
                        for i7 in range(top_words_limit[6][0], top_words_limit[6][1]+1):
                            for i8 in range(top_words_limit[7][0], top_words_limit[7][1]+1):
                                tmp_list = [1, 1] + [i1, i2, i3, i4, i5, i6, i7, i8]
                                new_code = part1 + our_data_head + '{' + ','.join([str(ti) for ti in tmp_list]) + '}' + our_data_tail + part2
                                open('keygen.c', 'w').write(new_code)
                                tmp = os.popen(compile_command1)
                                print('now list:', tmp_list)
                                for _ in range(128):
                                    tmp_out = tmp.readline()
                                    if 'win' in tmp_out:
                                        start_time = time.time()
                                        for _ in range(1000):
                                            tmp = os.popen(compile_command2)
                                        end_time = time.time()
                                        print(end_time - start_time)
                                        if(end_time-start_time < min_time):
                                            min_time = end_time-start_time
                                            print("win!!!!!!!!!!!!!!!!!")
                                            open('all_list.txt', 'a').write(str(tmp_list) + '\n')
                                            break
                                        else:
                                            break
                                    # if 'reduce err 2' in tmp_out:
                                    #     print("Err 22222!!!!!!!!!!")
                                    #     break
                                    # print(tmp_out)
                                      
