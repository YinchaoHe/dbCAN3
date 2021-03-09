import csv


def main():
    with open('prodigal.gff', 'r') as f:
        data = f.readlines()
    f.close()

    with open('overview.txt', 'r') as f:
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
            if new_id in ref_id:
                print(new_id)
                new_brief = brief.replace(columns[0], new_id, 1) + '\n'
                result.append(new_brief)

    with open('filtered_overv.gff', 'w') as file:
        file.writelines(result)
    file.close()

if __name__ == '__main__':
    main()