import numpy as np
import scipy as sp
import math


def hyperbolic_transform(distance):
    try:
        hyperbolic_distance = math.log((1 + distance) / (1 - distance))
    except ValueError:
        print(distance)
        print(1 + distance) / (1 - distance)
        raise ValueError

    return hyperbolic_distance


def hellinger(x, y):
    """Computes the Hellinger distance between two sum-normalized vectors

    :math:`d_h(x,y) = \\frac{1}{2} \sqrt{\sum_{i} \\left(\sqrt{x_i} - \sqrt{y_i}\\right)^2}`

    Equivalently, for sum-normalized vectors

    :math:`d_h(x,y) =  1 - \sum_{i} \\left(\sqrt{x_i} - \sqrt{y_i}\\right)`

    :param x: The first vector
    :param y: The second vector
    :type x: ndarray
    :type y: ndarray
    :returns: ``hellinger(x,y) = np.sqrt( 1 - np.dot(np.sqrt(x),np.sqrt(y)) )``
    """

    # hellinger_distance = sp.spatial.distance.euclidean(np.sqrt(x),np.sqrt(y))
    squared_hellinger_distance = 1 - np.dot(np.sqrt(x), np.sqrt(y))
    hellinger_distance = np.sqrt(squared_hellinger_distance)
    return hellinger_distance


def hellinger_hyp(x, y):
    """Computes the hyperbolic Hellinger distance metric between two vectors

    .. math::

       d_{hy} = \log \\left(\\frac{1 + d_h}{1 - d_h}\\right)

    :param x: The first vector
    :param y: The second vector
    :type x: ndarray
    :type y: ndarray
    :returns: ``hellinger_hyp(x,y) = math.log( (1+hellinger(x,y))/(1 - hellinger(x,y)) )``
    """
    hellinger_distance = hellinger(x, y)
    try:
        hellinger_hyperbolic_distance = hyperbolic_transform(hellinger_distance)
    except:
        print(hellinger_distance)
    return hellinger_hyperbolic_distance


def jeffries(x, y):
    """Computes the Jeffries divergence between two vectors

    .. math::

       d_{SKL}(x,y) = \\frac{1}{2} \\left( \sum_{i} \\left( x_i - y_i \\right) \log \\frac{x_i}{y_i}  \\right)


    :param x: The first vector
    :param y: The second vector
    :type x: ndarray
    :type y: ndarray
    :returns: ``jeffries(x,y) = (1/2) * ( entropy(x,y) + entropy(y,x) )``

    """
    entropy_xy = sp.stats.entropy(x, y)
    entropy_yx = sp.stats.entropy(y, x)
    jeffries_entropy = (entropy_xy + entropy_yx) / 2.0
    return jeffries_entropy


def jensen(x, y):
    """Computes the Jensen-Shannon metric between two vectors

    .. math::

       d^2_j = \sum_{i} \\left( x_i \log_2 \\frac{2x_i}{x_i + y_i} + y_i \log_2 \\frac{2y_i}{x_i + y_i} \\right)

    :param x: The first vector
    :param y: The second vector
    :type x: ndarray
    :type y: ndarray
    :returns: ``jensen(x,y) = entropy(x,m) + entropy(y,m)`` where ``m = (1/2) * (x + y)``

    """
    mean = (x + y) / 2.0
    entropy_xm = sp.stats.entropy(x, mean)
    entropy_ym = sp.stats.entropy(y, mean)
    jensen_entropy = (entropy_xm + entropy_ym) / (
        2.0 * math.log(2)
    )  # Log 2 is to go from nats to bits
    return jensen_entropy


def jensen_distance(x, y):
    """Computes the Jensen-Shannon distance between two vectors

    .. math::

       d_j = \sum_{i} \\left( x_i \log \\frac{2x_i}{x_i + y_i} + y_i \log \\frac{2y_i}{x_i + y_i} \\right)

    :param x: The first vector
    :param y: The second vector
    :type x: ndarray
    :type y: ndarray
    :returns: ``jensen_distance(x,y) = sqrt(entropy(x,m) + entropy(y,m))``

    where ``m = (1/2) * (x + y)``

    """
    jensen_entropy = jensen(x, y)
    jensen_distance = np.sqrt(jensen_entropy)
    return jensen_distance


def jensen_hyp(x, y):
    """Computes the hyperbolic Jensen-Shannon metric between two vectors

    .. math::

       d_{hy} = \log \\left(\\frac{1 + d_j}{1 - d_j}\\right)

    :param x: The first vector
    :param y: The second vector
    :type x: ndarray
    :type y: ndarray
    :returns: ``jensen_hyp(x,y) = math.log( (1+jensen_distance(x,y))/(1 - jensen_distance(x,y)) )``
    """
    jensen_entropy = jensen(x, y)
    jensen_distance = np.sqrt(jensen_entropy)
    jensen_hyperbolic_distance = hyperbolic_transform(jensen_distance)
    return jensen_hyperbolic_distance
