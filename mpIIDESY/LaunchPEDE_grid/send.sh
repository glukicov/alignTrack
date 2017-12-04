jobsub_submit -N 10000 -G gm2 -M file:///gm2/app/users/glukicov/submitGridJobs/grid_submit.sh --expected-lifetime=1h --memory=500MB --disk=2GB --resource-provides=usage_model=OFFSITE --role=Analysis

# jobsub_submit -N 1 -G gm2 -M file:///gm2/app/users/glukicov/submitGridJobs/grid_submit.sh --tar_file_name=///pnfs/gm2/scratch/users/glukicov/tracker.tar --expected-lifetime=1h \ --memory=500MB --disk=2GB \ --resource-provides=usage_model=OFFSITE \ --role=Analysis 
# jobsub_submit -N 1 -G gm2  -M file:///gm2/app/users/glukicov/submitGridJobs/grid_submit.sh file2:///pnfs/gm2/scratch/users/glukicov/tracker.tar --expected-lifetime=1h \ --memory=500MB --disk=2GB \ --resource-provides=usage_model=OFFSITE \ --role=Analysis 
# jobsub_submit -N 1 -G gm2  -M file:///gm2/app/users/glukicov/submitGridJobs/grid_submit.sh --f=///pnfs/gm2/scratch/users/glukicov/tracker.tar --expected-lifetime=1h \ --memory=500MB --disk=2GB \ --resource-provides=usage_model=OFFSITE \ --role=Analysis 
