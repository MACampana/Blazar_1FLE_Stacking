#!/bin/sh

#Arguments 
#$1 = Type of job (trials, erange-trials, sensdisc, erange-sens, diffsens, bias)
#$2 = hemisphere  (both, north, south)
#$3 = low or high cutoff for erange jobs (low, high, both)
#$4 = energy in GeV for erange jobs (number). Lower limit if $3=both.
#$5 = high energy limit for $3=both

mkdir '/scratch/mcampana/'

mkdir '/scratch/mcampana/outputs/'
mkdir '/scratch/mcampana/errors/'
mkdir '/scratch/mcampana/logs/'
mkdir '/scratch/mcampana/job_files/'

#rm -r "/scratch/mcampana/job_files/"*

mkdir '/scratch/mcampana/job_files/execs/'
mkdir "/scratch/mcampana/job_files/execs/$1/"
mkdir "/scratch/mcampana/job_files/execs/$1/$2"

mkdir '/scratch/mcampana/job_files/subs/'
mkdir "/scratch/mcampana/job_files/subs/$1/"
mkdir "/scratch/mcampana/job_files/subs/$1/$2"

mkdir '/scratch/mcampana/job_files/dags/'
mkdir "/scratch/mcampana/job_files/dags/$1/"
mkdir "/scratch/mcampana/job_files/dags/$1/$2"

gammas=(1.75 2.0 2.25 2.5 2.75 3.0 3.25)
weights=('equal' 'flux')

if [ $1 == 'trials' ]; then

    dag_path="/scratch/mcampana/job_files/dags/$1/$2/Dagman_$1_weights_gammas_$2_$3$4.dag"
    touch ${dag_path}
    
    COUNTER=({1..1000..1})
    for w in ${weights[@]}; do
        for g in ${gammas[@]}; do
            for c in ${COUNTER[@]}; do
                exec_path="/scratch/mcampana/job_files/execs/$1/$2/Do_$1_${w}weight_gamma${g}_seed${c}_$2.sh"
                touch ${exec_path}
                echo "#!/bin/sh" >> ${exec_path}
                echo "/data/user/mcampana/analysis/Blazar_1FLE/BlazarAnalysis.py -w ${w} -g ${g} -n 100 --cpus 1 -c --seed ${c} --hemi $2 --dataset ps-v3p2" >> ${exec_path}
    
                sub_path="/scratch/mcampana/job_files/subs/$1/$2/Submit_$1_${w}weight_gamma${g}_seed${c}_$2.submit"
                touch ${sub_path}
                echo "executable = ${exec_path}" >> ${sub_path}
                echo "output = /scratch/mcampana/outputs/$1_${w}weight_gamma${g}_seed${c}_$2.out" >> ${sub_path}
                echo "error = /scratch/mcampana/errors/$1_${w}weight_gamma${g}_seed${c}_$2.err" >> ${sub_path}
                echo "log = /scratch/mcampana/logs/$1_${w}weight_gamma${g}_seed${c}_$2.log" >> ${sub_path}        
                echo "getenv = true" >> ${sub_path}
                echo "universe = vanilla" >> ${sub_path}
                echo "notifications = never" >> ${sub_path}
                echo "should_transfer_files = YES" >> ${sub_path}
                echo "request_memory = 8000" >> ${sub_path}
                echo "queue 1" >> ${sub_path}
    
                echo "JOB $1.${w}.${g}.${c}.$2 ${sub_path}" >> ${dag_path}
    
            done
        done
    done
fi

if [ $1 == 'erange-trials' ]; then
    
    mkdir "/scratch/mcampana/job_files/execs/$1/$2/$3"
    mkdir "/scratch/mcampana/job_files/execs/$1/$2/$3/$4"
    
    mkdir "/scratch/mcampana/job_files/dags/$1/$2/$3"
    mkdir "/scratch/mcampana/job_files/dags/$1/$2/$3/$4"
    
    mkdir "/scratch/mcampana/job_files/subs/$1/$2/$3"
    mkdir "/scratch/mcampana/job_files/subs/$1/$2/$3/$4"
    
    dag_path="/scratch/mcampana/job_files/dags/$1/$2/$3/$4/Dagman_$1_weights_gammas_$2_$3$4.dag"
    touch ${dag_path}
    
    COUNTER=({1..100..1})
    
    g=2.0
    for w in ${weights[@]}; do
        for c in ${COUNTER[@]}; do
            exec_path="/scratch/mcampana/job_files/execs/$1/$2/$3/$4/Do_$1_${w}weight_$3$4_seed${c}_$2.sh"
            touch ${exec_path}
            echo "#!/bin/sh" >> ${exec_path}
            
            if [ $3 == 'low' ]; then
                echo "/data/user/mcampana/analysis/Blazar_1FLE/BlazarAnalysis.py -w ${w} -g ${g} -n 100 --cpus 1 -c --seed ${c} --hemi $2 --dataset ps-v3p2 --no-bg --ecut --elo $4" >> ${exec_path}
            fi
            if [ $3 == 'high' ]; then
                echo "/data/user/mcampana/analysis/Blazar_1FLE/BlazarAnalysis.py -w ${w} -g ${g} -n 100 --cpus 1 -c --seed ${c} --hemi $2 --dataset ps-v3p2 --no-bg --ecut --ehi $4" >> ${exec_path}
            fi
            if [ $3 == 'both' ]; then
                echo "/data/user/mcampana/analysis/Blazar_1FLE/BlazarAnalysis.py -w ${w} -g ${g} -n 100 --cpus 1 -c --seed ${c} --hemi $2 --dataset ps-v3p2 --no-bg --ecut --elo $4 --ehi $5" >> ${exec_path}
            fi
            
            sub_path="/scratch/mcampana/job_files/subs/$1/$2/$3/$4/Submit_$1_${w}weight_$3$4_seed${c}_$2.submit"
            touch ${sub_path}
            echo "executable = ${exec_path}" >> ${sub_path}
            echo "output = /scratch/mcampana/outputs/$1_${w}weight_$3$4_seed${c}_$2.out" >> ${sub_path}
            echo "error = /scratch/mcampana/errors/$1_${w}weight_$3$4_seed${c}_$2.err" >> ${sub_path}
            echo "log = /scratch/mcampana/logs/$1_${w}weight_$3$4_seed${c}_$2.log" >> ${sub_path}        
            echo "getenv = true" >> ${sub_path}
            echo "universe = vanilla" >> ${sub_path}
            echo "notifications = never" >> ${sub_path}
            echo "should_transfer_files = YES" >> ${sub_path}
            echo "request_memory = 8000" >> ${sub_path}
            echo "queue 1" >> ${sub_path}
    
            echo "JOB $1.${w}.${g}.${c}.$2.$3.$4 ${sub_path}" >> ${dag_path}
    
        done
    done
fi

if [ $1 == 'sensdisc' ]; then
    
    dag_path="/scratch/mcampana/job_files/dags/$1/$2/Dagman_$1_weights_gammas_$2_$3$4.dag"
    touch ${dag_path}

    exec_path="/scratch/mcampana/job_files/execs/$1/$2/Do_$1_$2.sh"
    touch ${exec_path}
    echo "#!/bin/sh" >> ${exec_path}
    echo "/data/user/mcampana/analysis/Blazar_1FLE/BlazarAnalysis.py -w ${weights[@]} -g ${gammas[@]} -l -s -d --hemi $2 --dataset ps-v3p2" >> ${exec_path}

    sub_path="/scratch/mcampana/job_files/subs/$1/$2/Submit_$1_$2.submit"
    touch ${sub_path}
    echo "executable = ${exec_path}" >> ${sub_path}
    echo "output = /scratch/mcampana/outputs/$1_$2.out" >> ${sub_path}
    echo "error = /scratch/mcampana/errors/$1_$2.err" >> ${sub_path}
    echo "log = /scratch/mcampana/logs/$1_$2.log" >> ${sub_path}        
    echo "getenv = true" >> ${sub_path}
    echo "universe = vanilla" >> ${sub_path}
    echo "notifications = never" >> ${sub_path}
    echo "should_transfer_files = YES" >> ${sub_path}
    echo "request_memory = 10000" >> ${sub_path}
    echo "queue 1" >> ${sub_path}

    echo "JOB $1.$2 ${sub_path}" >> ${dag_path}  
        
fi

if [ $1 == 'erange-sens' ]; then

    mkdir "/scratch/mcampana/job_files/execs/$1/$2/$3"
    mkdir "/scratch/mcampana/job_files/execs/$1/$2/$3/$4"
    
    mkdir "/scratch/mcampana/job_files/dags/$1/$2/$3"
    mkdir "/scratch/mcampana/job_files/dags/$1/$2/$3/$4"
    
    mkdir "/scratch/mcampana/job_files/subs/$1/$2/$3"
    mkdir "/scratch/mcampana/job_files/subs/$1/$2/$3/$4"
    
    dag_path="/scratch/mcampana/job_files/dags/$1/$2/$3/$4/Dagman_$1_weights_gammas_$2_$3$4.dag"
    touch ${dag_path}
    
    exec_path="/scratch/mcampana/job_files/execs/$1/$2/$3/$4/Do_$1_$3$4_$2.sh"
    touch ${exec_path}
    echo "#!/bin/sh" >> ${exec_path}
    
    if [ $3 == 'low' ]; then
        echo "/data/user/mcampana/analysis/Blazar_1FLE/BlazarAnalysis.py -w ${weights[@]} -g 2.0 -l -s -d --hemi $2 --dataset ps-v3p2 --ecut --elo $4" >> ${exec_path}
    fi
    if [ $3 == 'high' ]; then
        echo "/data/user/mcampana/analysis/Blazar_1FLE/BlazarAnalysis.py -w ${weights[@]} -g 2.0 -l -s -d --hemi $2 --dataset ps-v3p2 --ecut --ehi $4" >> ${exec_path}
    fi
    if [ $3 == 'both' ]; then
        echo "/data/user/mcampana/analysis/Blazar_1FLE/BlazarAnalysis.py -w ${weights[@]} -g 2.0 -l -s -d --hemi $2 --dataset ps-v3p2 --ecut --elo $4 --ehi $5" >> ${exec_path}
    fi

    sub_path="/scratch/mcampana/job_files/subs/$1/$2/$3/$4/Submit_$1_$3$4_$2.submit"
    touch ${sub_path}
    echo "executable = ${exec_path}" >> ${sub_path}
    echo "output = /scratch/mcampana/outputs/$1_$3$4_$2.out" >> ${sub_path}
    echo "error = /scratch/mcampana/errors/$1_$3$4_$2.err" >> ${sub_path}
    echo "log = /scratch/mcampana/logs/$1_$3$4_$2.log" >> ${sub_path}        
    echo "getenv = true" >> ${sub_path}
    echo "universe = vanilla" >> ${sub_path}
    echo "notifications = never" >> ${sub_path}
    echo "should_transfer_files = YES" >> ${sub_path}
    echo "request_memory = 10000" >> ${sub_path}
    echo "queue 1" >> ${sub_path}

    echo "JOB $1.$3.$4.$2 ${sub_path}" >> ${dag_path}  
        
fi

if [ $1 == 'diffsens' ]; then
    
    dag_path="/scratch/mcampana/job_files/dags/$1/$2/Dagman_$1_weights_gammas_$2_$3$4.dag"
    touch ${dag_path}
    
    exec_path="/scratch/mcampana/job_files/execs/$1/$2/Do_$1_$2.sh"
    touch ${exec_path}
    echo "#!/bin/sh" >> ${exec_path}
    echo "/data/user/mcampana/analysis/Blazar_1FLE/BlazarAnalysis.py -w ${weights[@]} -g 2.0 -l -s --diff-sens --hemi $2 --dataset ps-v3p2" >> ${exec_path}

    sub_path="/scratch/mcampana/job_files/subs/$1/$2/Submit_$1_$2.submit"
    touch ${sub_path}
    echo "executable = ${exec_path}" >> ${sub_path}
    echo "output = /scratch/mcampana/outputs/$1_$2.out" >> ${sub_path}
    echo "error = /scratch/mcampana/errors/$1_$2.err" >> ${sub_path}
    echo "log = /scratch/mcampana/logs/$1_$2.log" >> ${sub_path}        
    echo "getenv = true" >> ${sub_path}
    echo "universe = vanilla" >> ${sub_path}
    echo "notifications = never" >> ${sub_path}
    echo "should_transfer_files = YES" >> ${sub_path}
    echo "request_memory = 10000" >> ${sub_path}
    echo "queue 1" >> ${sub_path}

    echo "JOB $1.$2 ${sub_path}" >> ${dag_path}  
        
fi

if [ $1 == 'bias' ]; then

    gammas=(2.0 2.5 3.0)

    dag_path="/scratch/mcampana/job_files/dags/$1/$2/Dagman_$1_weights_gammas_$2_$3$4.dag"
    touch ${dag_path}
    
    for w in ${weights[@]}; do
        for g in ${gammas[@]}; do
            exec_path="/scratch/mcampana/job_files/execs/$1/$2/Do_$1_${w}_${g}_$2.sh"
            touch ${exec_path}
            echo "#!/bin/sh" >> ${exec_path}
            echo "/data/user/mcampana/analysis/Blazar_1FLE/BlazarAnalysis.py -w ${w} -g ${g} -l -b --hemi $2 --no-poisson --dataset ps-v3p2" >> ${exec_path}

            sub_path="/scratch/mcampana/job_files/subs/$1/$2/Submit_$1_${w}_${g}_$2.submit"
            touch ${sub_path}
            echo "executable = ${exec_path}" >> ${sub_path}
            echo "output = /scratch/mcampana/outputs/$1_${w}_${g}_$2.out" >> ${sub_path}
            echo "error = /scratch/mcampana/errors/$1_${w}_${g}_$2.err" >> ${sub_path}
            echo "log = /scratch/mcampana/logs/$1_${w}_${g}_$2.log" >> ${sub_path}        
            echo "getenv = true" >> ${sub_path}
            echo "universe = vanilla" >> ${sub_path}
            echo "notifications = never" >> ${sub_path}
            echo "should_transfer_files = YES" >> ${sub_path}
            echo "request_memory = 8000" >> ${sub_path}
            echo "queue 1" >> ${sub_path}

            echo "JOB $1.${w}.${g}.$2 ${sub_path}" >> ${dag_path}  
        
        done
    done        
fi

runThis="/scratch/mcampana/job_files/SubmitMyJobs_$1_$2_$3$4.sh"
touch ${runThis}
echo "#!/bin/sh" >> ${runThis}
echo "condor_submit_dag -maxjobs 500 ${dag_path}" >> ${runThis}

