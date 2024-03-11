set -ex

cd ../../lavaan_fitmodel
tar -zcvf lavaan.tar.gz lavaan
cd ../TRAINING/redux10point2

zone=us-central1-a
prefix=ccarey
script_name=fitmeasures

name=$prefix-$script_name-training
gcloud compute instances create --zone $zone --machine-type m1-ultramem-80 --boot-disk-size=500GB --image cfa-benchmarking $name &

wait

sleep 5 # give the instances a moment to fully come alive

gcloud compute scp --zone $zone $script_name.r ../../lavaan_fitmodel/lavaan.tar.gz fit_full.Rds $name:. 
