import sys

ikbkg_true_start = 154547620
ikbkg_match_p1_start = 154556163
ikbkg_match_p1_end = 154565046

# notice how they are flipped
ikbkgp1_end = 154639391
ikbkgp1_start = 154648275

ikbkgp1_len = 8885

def get_matching_ikbkgp1_pos(position):
    global ikbkg_match_p1_start
    global ikbkgp1_start

    # there is an gap in the reference just after this position
    # offset 15544
    if position <= 154563013:
        offset = position - ikbkg_match_p1_start
    else:
        offset = position - ikbkg_match_p1_start + 1

    #print(offset)
    return ikbkgp1_start - offset

def get_matching_ikbkg_pos(position):
    global ikbkg_match_p1_start
    global ikbkgp1_start

    if position >= 154641425:
        offset = ikbkgp1_start - position
    elif position == 154641424:
        raise Exception("There is an extra \"A\" at this exact IKBKGP1 position. Case not handled yet.")
    else:
        offset = ikbkgp1_start - position - 1 # 154641423 (IKBKGP1) == 154563014 (IKBKG)

    return ikbkg_match_p1_start + offset

print(get_matching_ikbkgp1_pos(154560403))
