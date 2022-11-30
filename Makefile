unit:
	brawn --help > /dev/null
	pytest

lint:
	mypy brawn
	pylint -E brawn
	pylint -rn brawn

coverage:
	coverage run --source . -m pytest
	coverage report
	coverage html
