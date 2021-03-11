def gff_filter(gff_file, overview_file, output_file):
    with open(gff_file, 'r') as f:
        data = f.readlines()
    f.close()

    with open(overview_file, 'r') as f:
        raw_refs = f.readlines()
    f.close()


    ref_id =[]
    for raw_ref in raw_refs[1:]:
        id = raw_ref.split()[0]
        ref_id.append(id)


    result = []
    for row in data:
        if row.startswith('#'):
           continue
        else:
            brief = row.split('\n')[0]
            columns = brief.split()
            des = columns[-1]
            id = des.split(';partial')[0]
            id_suffix = id.split('_')[1]
            new_id = columns[0] + '_' + id_suffix
            output_id = 'ID='+new_id
            if new_id in ref_id:
                print(new_id)
                new_brief = brief.replace(id, output_id, 1) + '\n'
                result.append(new_brief)

    with open(output_file, 'w') as file:
        file.writelines(result)
    file.close()

def check_result(new_id):
    with open('uniInput', 'r') as f:
        raw_checks = f.readlines()
    f.close()

    checks = []
    for raw_check in raw_checks:
        if raw_check.startswith('>'):
            raw_id = raw_check.split(' # ')[0]
            id = raw_id.split('>')[1]
            checks.append(id)

    if new_id in checks:
        print('checked:' + new_id)
    else:
        print(new_id)

if __name__ == '__main__':
    gff_filter(gff_file='prodigal.gff', overview_file='overview.txt', output_file='filtered_overv.gff')