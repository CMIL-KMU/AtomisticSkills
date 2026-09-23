from src.utils.paper_downloader import _is_url_blocked


def test_is_url_blocked_filtering():
    assert (
        _is_url_blocked("https://link.aps.org/pdf/10.1103/PhysRevLett.120.145301")
        is True
    )
    assert _is_url_blocked("https://pubs.acs.org/doi/pdf/10.1021/jacs.3c01234") is True
    assert (
        _is_url_blocked(
            "https://www.sciencedirect.com/science/article/pii/S092702562100520X"
        )
        is True
    )
    assert (
        _is_url_blocked(
            "https://onlinelibrary.wiley.com/doi/pdf/10.1002/anie.202300001"
        )
        is True
    )

    # Safe repositories must NOT be blocked
    assert _is_url_blocked("https://arxiv.org/pdf/2301.00001.pdf") is False
    assert (
        _is_url_blocked(
            "https://europepmc.org/backend/ptpmcrender.fcgi?accid=PMC4336040&blobtype=pdf"
        )
        is False
    )
    assert (
        _is_url_blocked(
            "https://dspace.mit.edu/bitstreams/b6249448-f14f-4805-bf14-3e013cd01173/download"
        )
        is False
    )
    assert (
        _is_url_blocked("https://zenodo.org/record/123456/files/article.pdf") is False
    )
    assert _is_url_blocked("https://www.osti.gov/servlets/purl/1524040") is False
