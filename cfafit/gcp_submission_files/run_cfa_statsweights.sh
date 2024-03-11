set -ex

cd ../../lavaan_statsweights
tar -zcvf lavaan.tar.gz lavaan
cd ../TRAINING/redux10point2

zone=us-central1-a
prefix=ccarey
script_name=cfa2
cores=80

name=$prefix-$script_name-training
gcloud compute instances create --zone $zone --machine-type m1-ultramem-$cores --boot-disk-size=500GB --image cfa-benchmarking $name &

wait

sleep 5 # give the instances a moment to fully come alive

gcloud compute scp --zone $zone $script_name.r ../../lavaan_statsweights/lavaan.tar.gz ../../imp1.Rds ../../full.txt ../train_withdraw_outliers.txt $name:.
