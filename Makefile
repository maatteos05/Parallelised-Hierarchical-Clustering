.PHONY: all average single clean

all: average single

average:
	$(MAKE) -C average-link

single:
	$(MAKE) -C single-link

clean:
	$(MAKE) -C average-link clean
	$(MAKE) -C single-link clean
