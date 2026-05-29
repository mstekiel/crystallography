# noxfile.py
import nox

PYTHON_VERSIONS = ["3.10", "3.11", "3.12"]
NUMPY_VERSIONS = ["1.23", "1.26", "2.0", "2.4"]

@nox.session
def tests(session):
    """Run pytest with coverage."""
    session.install("pytest", "pytest-cov")
    session.install("-e", ".")
    session.run("pytest")

@nox.session
def coverage(session):
    """Generate HTML coverage report."""
    session.install("pytest", "pytest-cov")
    session.install("-e", ".")
    session.run("pytest", "--cov=bragg", "--cov-report=html")

@nox.session(python=PYTHON_VERSIONS)
@nox.parametrize("numpy", NUMPY_VERSIONS)
def version_matrix(session, numpy):
    """Test across Python and numpy versions."""
    session.install(f"numpy=={numpy}")
    session.install("pytest")
    session.install("-e", ".")
    session.run("pytest")