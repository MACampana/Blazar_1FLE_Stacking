#!/bin/sh
mkdir '/scratch/mcampana/'

mkdir '/scratch/mcampana/outputs/'
mkdir '/scratch/mcampana/errors/'
mkdir '/scratch/mcampana/logs/'
mkdir '/scratch/mcampana/job_files/'

rm -r "/scratch/mcampana/job_files/"*

mkdir '/scratch/mcampana/job_files/execs/'
mkdir '/scratch/mcampana/job_files/subs/'
mkdir '/scratch/mcampana/job_files/dags/'

gammas=(1.75 2.0 2.25 2.5 2.75 3.0 3.25)
weights=('equal' 'flux')

dag_path="/scratch/mcampana/job_files/dags/Dagman_$1_weights_gammas.dag"
touch ${dag_path}


if [ $1 == 'trials' ]; then
    COUNTER=({1..100..1})
    for w in ${weights[@]}; do
        for g in ${gammas[@]}; do
            for c in ${COUNTER[@]}; do
                exec_path="/scratch/mcampana/job_files/execs/Do_$1_${w}weight_gamma${g}_seed${c}.sh"
                touch ${exec_path}
                echo "#!/bin/sh" >> ${exec_path}
                echo "/data/user/mcampana/analysis/Blazar_1FLE/BlazarAnalysis.py -w ${w} -g ${g} -n 100 --cpus 1 -c --seed ${c}" >> ${exec_path}
    
                sub_path="/scratch/mcampana/job_files/subs/Submit_$1_${w}weight_gamma${g}_seed${c}.submit"
                touch ${sub_path}
                echo "executable = ${exec_path}" >> ${sub_path}
                echo "output = /scratch/mcampana/outputs/$1_${w}weight_gamma${g}_seed${c}.out" >> ${sub_path}
                echo "error = /scratch/mcampana/errors/$1_${w}weight_gamma${g}_seed${c}.err" >> ${sub_path}
                echo "log = /scratch/mcampana/logs/$1_${w}weight_gamma${g}_seed${c}.log" >> ${sub_path}        
                echo "getenv = true" >> ${sub_path}
                echo "universe = vanilla" >> ${sub_path}
                echo "notifications = never" >> ${sub_path}
                echo "should_transfer_files = YES" >> ${sub_path}
                echo "request_memory = 8000" >> ${sub_path}
                echo "queue 1" >> ${sub_path}
    
                echo "JOB $1.${w}.${g}.${c} ${sub_path}" >> ${dag_path}
    
            done
        done
    done
fi


if [ $1 == 'sensdisc' ]; then
    exec_path="/scratch/mcampana/job_files/execs/Do_$1.sh"
    touch ${exec_path}
    echo "#!/bin/sh" >> ${exec_path}
    echo "/data/user/mcampana/analysis/Blazar_1FLE/BlazarAnalysis.py -w ${weights[@]} -g ${gammas[@]} -l -s -i " >> ${exec_path}

    sub_path="/scratch/mcampana/job_files/subs/Submit_$1_.submit"
    touch ${sub_path}
    echo "executable = ${exec_path}" >> ${sub_path}
    echo "output = /scratch/mcampana/outputs/$1_.out" >> ${sub_path}
    echo "error = /scratch/mcampana/errors/$1_.err" >> ${sub_path}
    echo "log = /scratch/mcampana/logs/$1_.log" >> ${sub_path}        
    echo "getenv = true" >> ${sub_path}
    echo "universe = vanilla" >> ${sub_path}
    echo "notifications = never" >> ${sub_path}
    echo "should_transfer_files = YES" >> ${sub_path}
    echo "request_memory = 10000" >> ${sub_path}
    echo "queue 1" >> ${sub_path}

    echo "JOB $1 ${sub_path}" >> ${dag_path}  
        
fi

if [ $1 == 'bias' ]; then
    exec_path="/scratch/mcampana/job_files/execs/Do_$1.sh"
    touch ${exec_path}
    echo "#!/bin/sh" >> ${exec_path}
    echo "/data/user/mcampana/analysis/Blazar_1FLE/BlazarAnalysis.py -w ${weights[@]} -g ${gammas[@]} -b " >> ${exec_path}

    sub_path="/scratch/mcampana/job_files/subs/Submit_$1_.submit"
    touch ${sub_path}
    echo "executable = ${exec_path}" >> ${sub_path}
    echo "output = /scratch/mcampana/outputs/$1_.out" >> ${sub_path}
    echo "error = /scratch/mcampana/errors/$1_.err" >> ${sub_path}
    echo "log = /scratch/mcampana/logs/$1_.log" >> ${sub_path}        
    echo "getenv = true" >> ${sub_path}
    echo "universe = vanilla" >> ${sub_path}
    echo "notifications = never" >> ${sub_path}
    echo "should_transfer_files = YES" >> ${sub_path}
    echo "request_memory = 8000" >> ${sub_path}
    echo "queue 1" >> ${sub_path}

    echo "JOB $1 ${sub_path}" >> ${dag_path}  
        
fi

runThis="/scratch/mcampana/job_files/SubmitMyJobs.sh"
touch ${runThis}
echo "#!/bin/sh" >> ${runThis}
echo "condor_submit_dag ${dag_path}" >> ${runThis}

