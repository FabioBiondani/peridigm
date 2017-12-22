for i in `find src -name "*"`; do \
    if [ -f $i ]; then \
        sed s/std::tr1/std/g $i > tmp; \
        mv tmp $i; \
    fi; \
done
