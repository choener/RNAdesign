
I=1

for i in $1/*constraints
do
  cp $i $2/$I.constraints
  I=$(( $I+1 ))
done

