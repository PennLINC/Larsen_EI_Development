import flywheel
import os

fw = flywheel.Client()

project = fw.lookup('bbl/PNC_CS_810336')
sessions = project.sessions()
subjects = project.subjects()

analysis_str = 'IDEMO_ACOMPCOR_GSR'#This is the newest version. Previously we used 'acompcor' and 'SDK_TASK'

for sub in subjects:
    """Loop over subjects and get each session"""
    sub_label = sub.label
    # if not "14665" in sub_label:
#         continue
    
    for ses in sub.sessions(): 
        ses_label = ses.label
        """Get the analyses for that session"""
        full_ses = fw.get(ses.id)
        these_analyses = [ana for ana in full_ses.analyses if analysis_str in ana.label]
        these_analyses_labs = [ana.label for ana in full_ses.analyses if analysis_str in ana.label] 
        if len(these_analyses)<1:
             print('No analyses {} {}'.format(sub_label,ses_label))
             continue
        for this_ana in these_analyses:
            if not this_ana.files:
                continue

            outputs = [f for f in this_ana.files if f.name.endswith('.zip')
                and not f.name.endswith('.html.zip')]
            output = outputs[0]
            
            ana_label = this_ana.label.split(' ')[0]   
            
            dest = '/cbica/projects/alpraz_EI/data/PNC/{}/{}/{}/'.format(analysis_str,sub_label,ses_label)
            try:
                os.makedirs(dest)
            except OSError:
                print(dest+" exists")
            else: print("creating "+dest)
            dest_file = dest+output.name
            if not os.path.exists(dest_file):
                print("Downloading", dest_file)
                output.download(dest_file)
                print('Done')
