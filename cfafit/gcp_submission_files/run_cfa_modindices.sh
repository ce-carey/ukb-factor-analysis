set -ex

cd ../../lavaan_postfit_mod
tar -zcvf lavaan.tar.gz lavaan
cd ../TRAINING/redux10point2

zone=us-east1-b
prefix=ccarey
script_name=modindices

name=$prefix-$script_name-training
gcloud compute instances create --zone $zone --machine-type m1-ultramem-80 --boot-disk-size=500GB --image cfa-benchmarking $name &

wait

sleep 5 # give the instances a moment to fully come alive

gcloud compute scp --zone $zone $script_name.r ../../lavaan_postfit_mod/lavaan.tar.gz fit.Rds muthen1984.Rds lav_model_estimate.Rds $name:. 
