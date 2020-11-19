import os, sys

if sys.argv[1] == "-h":
    print("Below is a mockup of how to run the script:")
    print(
        "python3	b_cigar.py 	filename	 lower_conditional	 upper_conditional 	> 	output"
    )
    sys.exit()

filename = sys.argv[1]
lower_conditional = int(sys.argv[2])
upper_conditional = int(sys.argv[3])


for i, line in enumerate(open(filename)):
    if line.startswith("Qname"):
        print(line[:-1])
    else:
        l = line.strip().split("\t")
        if len(l) == 38 or len(l) == 40:

            b_cigar = l[16]

            if b_cigar == "NA":
                continue
            integer_list = []
            code_list = []
            temp_number = ""
            for c in b_cigar:
                if c.isdigit():
                    temp_number += c
                else:
                    integer_list.append(int(temp_number))
                    code_list.append(c)
                    temp_number = ""

            location_value = 0

            if sum(integer_list) <= upper_conditional:
                print(line[:-1])
                continue

            check = ""
            for i in range(len(integer_list)):
                if location_value + integer_list[i] < lower_conditional:
                    location_value += integer_list[i]
                elif (
                    location_value + integer_list[i] >= lower_conditional
                    and location_value + integer_list[i] <= upper_conditional
                ):
                    check += code_list[i]
                    location_value += integer_list[i]
                else:
                    # THIS IS WHERE IT SEEMS TO BE BROKEN
                    if not check and (code_list[i] == "M" or code_list[i] == "X"):
                        # print(line[:-1])
                        break
                    if "X" not in check:
                        print(line[:-1])
                        break
                    elif "X" in check:
                        if "D" in check or "I" in check:
                            print(line[:-1])
                            break
                        else:
                            break

        else:
            print(line[:-1])
