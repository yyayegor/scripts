"""Module for parsing results of NBO - quantum-chemical calculations
and getting info about second order perturbation energy for bonding
and unbonding interactions for bonds in molecule"""
import os
import re
from glob import glob


def lpline(analyzed_line, analyzed_bond):
    """Analysis of influence of concidered interaction on current
    bond, if there is a lone electron pair involved in interaction."""
    lpp = r'LP[^A-Za-z]*(\w+)\s+(\d+).+BD\*\( (\d+)\) (\w+)\s+(\d+)- (\w+)\s+(\d+)\s+([0-9.]+)'
    reg_line = re.search(lpp, analyzed_line)
    lpa = ''.join(reg_line.group(1, 2))
    acctp = reg_line.group(3)
    bd1 = ''.join(reg_line.group(4, 5))
    bd2 = ''.join(reg_line.group(6, 7))
    en = reg_line.group(8)
    if (lpa in analyzed_bond) and ((bd1 in analyzed_bond) or (bd2 in analyzed_bond)):
        stz = 'stabilizing'
    elif (bd1 in analyzed_bond) and (bd2 in analyzed_bond):
        stz = 'unbonding'
    else:
        stz = None
    if stz:
        return ('LP %s -> BD%s %s-%s : %s, %s for %s\n' %
                (lpa, acctp, bd1, bd2, en, stz, '-'.join(analyzed_bond)))
    else:
        return ''


def bdline(analyzed_line, analyzed_bond):
    """Analysis of influence of selected interaction on current
        bond, if there is no lone pairs."""
    bdp = r'BD \( (\d)\)\s(\w+)\s*(\d+)- (\w+)\s+(\d+).+' \
          r'BD\*\( (\d)\) (\w+)\s*(\d+)- (\w+)\s*(\d*)\s*([0-9.]*)'
    reg_line = re.search(bdp, analyzed_line)
    tp1 = reg_line.group(1)
    bd1 = ''.join(reg_line.group(2, 3))
    bd2 = ''.join(reg_line.group(4, 5))
    tp2 = reg_line.group(6)
    bd3 = ''.join(reg_line.group(7, 8))
    bd4 = ''.join(reg_line.group(9, 10))
    en = reg_line.group(11)
    sb = set(analyzed_bond)
    if ({bd1, bd2} == sb) or ({bd3, bd4} == sb):
        stz = 'unbonding'
    elif ({bd1, bd3} == sb) or ({bd1, bd4} == sb) \
            or ({bd2, bd3} == sb) or ({bd2, bd4} == sb):
        stz = 'stabilizing'
    else:
        stz = None
    if stz:
        return ('BD%s %s-%s -> BD%s %s-%s : %s, %s for %s\n' %
                (tp1, bd1, bd2, tp2, bd3, bd4, en, stz, '-'.join(analyzed_bond)))
    else:
        return ''


def run():
    filelist = [flnm for flnm in glob('*.nbo') if os.path.isfile(flnm)]
    if not filelist:
        input('No .nbo files found')
        exit()

    for filename in filelist:
        flnm = filename[:-3]
        with open(filename) as nbo:
            with open(flnm + 'sop', 'wt') as sop:
                flag1 = False
                flag2 = False
                for i in nbo:  # extracting second order perturbation data
                    if 'SECOND' in i:  # and writing it in .sop file
                        flag1 = True
                        continue
                    elif 'within' in i and flag1:
                        flag2 = True
                        continue
                    if flag1 and flag2 and 'Summary' in i:
                        break
                    if flag1 and flag2:
                        sop.write(i)

        with open(flnm + 'sop') as sop:
            atp = r'\s+\d+\.\s+BD \( \d\)\s+(\w+)\s+(\d+)- (\w+)\s+(\d+)'
            atoms = []
            bonds = []
            for strng in sop:
                g = re.match(atp, strng)
                if g:
                    a1, a2 = ''.join(g.group(1, 2)), ''.join(g.group(3, 4))
                    bond = (a1, a2)
                    if a1 not in atoms:
                        atoms.append(a1)
                    if a2 not in atoms:
                        atoms.append(a2)
                    if (bond not in bonds) and (bond[::-1] not in bonds):
                        bonds.append(bond)
        atoms.sort()
        n2o = {}
        if str(input('Would you like rename atoms? (Y/N): ')).lower() == 'y':
            rename = True
            for aa in atoms:
                n2o[aa] = str(input('%s : ' % aa)).capitalize()
        else:
            rename = False
            for aa in atoms:
                n2o[aa] = aa
        otpb = []
        for bond in bonds:
            with open(flnm + 'sop') as sop:
                for line in sop:
                    if ('RY' in line) or ('CR' in line) or line.isspace():
                        continue
                    if 'LP' in line:
                        otpb.append(lpline(line, bond))
                    else:
                        otpb.append(bdline(line, bond))
                otpb.append('\n')

        with open(flnm + 'otp', 'wt') as otp:
            if rename:
                for line in otpb:
                    for key in n2o:
                        line = line.replace(key, n2o[key])
                    otp.write(line)
            else:
                otp.writelines(otpb)


if __name__ == "__main__":
    run()
    input('Job complete')
